/* Exercise the real browser functions with an instrumented DOM. No screenshot
 * or physical validation is implied. Run with Node's standard library only. */
const assert = require('node:assert/strict');
const crypto = require('node:crypto');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const root = path.resolve(__dirname, '..');
const sourcePath = path.join(root, 'src/rocket_propulsion/web/app.js');
const source = fs.readFileSync(sourcePath, 'utf8');
class Element {
  constructor(tag) { this.tag=tag; this.children=[]; this.attributes={}; this.textContent=''; }
  setAttribute(key,value) { this.attributes[key]=String(value); }
  replaceChildren(...children) { this.children=children; }
  append(...children) { this.children.push(...children); }
}
const elements = new Map();
const document = {
  createElement: tag => new Element(tag),
  createElementNS: (namespace,tag) => new Element(tag),
  getElementById: id => {
    if (!elements.has(id)) elements.set(id,new Element('div'));
    return elements.get(id);
  },
};
const context = vm.createContext({document, console});
const prefix = source.slice(0,source.indexOf('function renderPolytropicCharts('));
assert.ok(prefix.includes('function renderLineChart('));
vm.runInContext(prefix, context);
const renderSource = source.slice(source.indexOf('function renderUncertainty('),source.indexOf('function analysisWorkspace('));
assert.ok(renderSource.startsWith('function renderUncertainty('));
vm.runInContext('let latestUncertaintyResult=null;\n'+renderSource,context);
const cycleStart=source.includes('function updateCycleControls(')
  ? source.indexOf('function updateCycleControls(') : source.indexOf('function renderEngineCycle(');
vm.runInContext('let latestCycleResult=null;\n'+source.slice(cycleStart,source.indexOf('function renderUncertainty(')),context);
const rows=[];
function check(name, action) {
  try { action(); rows.push({name,status:'passed'}); }
  catch(error) { rows.push({name,status:'failed',error:String(error)}); }
}
function chart(points) {
  const svg = new Element('svg');
  context.renderLineChart(svg,[{data:points,x:'x',y:'y'}],{xLabel:'x',yLabel:'y'});
  for (const child of svg.children) {
    for (const value of Object.values(child.attributes)) assert.doesNotMatch(value,/NaN|Infinity/);
    assert.doesNotMatch(child.textContent,/NaN|Infinity/);
  }
  const paths=svg.children.filter(child=>child.tag==='path');
  assert.ok(paths.length);
  return {svg,paths};
}
check('opposite finite extremes have finite coordinates and ordered endpoints',()=>{
  const {paths}=chart([{x:-Number.MAX_VALUE,y:-Number.MAX_VALUE},{x:Number.MAX_VALUE,y:Number.MAX_VALUE}]);
  const values=paths[0].attributes.d.match(/-?\d+\.\d+/g).map(Number);
  assert.ok(values[0]<values[2]); assert.ok(values[1]>values[3]);
});
check('constant maximum stays finite',()=>chart([{x:1,y:Number.MAX_VALUE},{x:2,y:Number.MAX_VALUE}]));
check('tiny nonconstant data occupies a visible range',()=>{
  const {paths}=chart([{x:1,y:1e-200},{x:2,y:2e-200}]);
  const values=paths[0].attributes.d.match(/-?\d+\.\d+/g).map(Number);
  assert.ok(Math.abs(values[1]-values[3])>100);
});
check('subnormal data and all-zero data stay finite',()=>{
  chart([{x:0,y:0},{x:Number.MIN_VALUE,y:Number.MIN_VALUE}]);
  chart([{x:0,y:0},{x:0,y:0}]);
});
check('invalid points break the line without malformed SVG',()=>{
  const {paths}=chart([{x:0,y:0},{x:1,y:NaN},{x:2,y:2}]);
  assert.equal((paths[0].attributes.d.match(/M/g)||[]).length,2);
});
check('declared maximum sample count renders without an argument-stack overflow',()=>{
  chart(Array.from({length:100000},(_,i)=>({x:i,y:i%100})));
});
check('null Pearson is explicitly undefined and parameter names are text',()=>{
  context.renderUncertainty({method:'latin-hypercube',output_percentiles:{p05:0,p95:1},
    sensitivity:[{parameter:'<img onerror=bad>',correlation:null,rank:null},
                 {parameter:'x',correlation:0,rank:1}],samples:[{output:0},{output:1}]});
  const list=document.getElementById('uncertainty-sensitivity');
  assert.equal(list.children.length,2);
  assert.match(list.children[0].textContent,/tanımsız.*sıralama yok/);
  assert.match(list.children[0].textContent,/<img onerror=bad>/);
  assert.equal(list.children[0].children.length,0);
  assert.match(list.children[1].textContent,/1\. x: Pearson r = 0/);
  assert.match(document.getElementById('uncertainty-inference').textContent,/güven aralığı değildir/);
});
check('pressure-fed shortage is visible and is never labelled achieved delivery',()=>{
  context.renderEngineCycle({cycle:'pressure-fed',pressure_fed_supply_sufficient:false,
    warnings:[{code:'insufficient_pressure_fed_supply',message:'insufficient'}]});
  assert.match(document.getElementById('cycle-scope-note').textContent,/karşılamıyor/);
});
check('external generator feed and untrusted warning text remain visible literal text',()=>{
  context.renderEngineCycle({cycle:'gas-generator',pressure_fed_supply_sufficient:null,
    inputs:{pump_feed_policy:'external-gg-feed'},warnings:[{code:'custom',message:'<img onerror=bad>'}]});
  const note=document.getElementById('cycle-scope-note');
  assert.match(note.textContent,/ayrı/);
  assert.match(note.textContent,/<img onerror=bad>/);
  assert.equal(note.children.length,0);
});
check('cycle selection disables inactive component fields and preserves their values',()=>{
  assert.equal(typeof context.updateCycleControls,'function');
  const names=['propellant_density_kg_m3','pump_efficiency','turbine_efficiency',
    'turbine_inlet_temperature_k','turbine_specific_heat_j_kg_k','turbine_gamma','turbine_pressure_ratio'];
  const fields=Object.fromEntries(names.map(name=>[name,{value:'42',disabled:false}]));
  fields.cycle={value:'pressure-fed'};fields.pump_feed_policy={disabled:false,value:'shared-pump'};
  const form={elements:fields};
  context.updateCycleControls(form);
  assert.ok(names.every(name=>fields[name].disabled));assert.ok(fields.pump_feed_policy.disabled);
  fields.cycle.value='gas-generator';context.updateCycleControls(form);
  assert.ok(names.every(name=>!fields[name].disabled&&fields[name].value==='42'));
  assert.equal(fields.pump_feed_policy.disabled,false);
  fields.cycle.value='staged-combustion';context.updateCycleControls(form);
  assert.ok(names.every(name=>!fields[name].disabled));assert.ok(fields.pump_feed_policy.disabled);
});
const report={scope:'Actual chart/render/control function executions with instrumented DOM; no browser layout or physics claim',
  source_sha256:crypto.createHash('sha256').update(source).digest('hex'),node:process.version,rows,
  passed:rows.every(row=>row.status==='passed')};
const target=process.argv[2];
if (target) fs.writeFileSync(target,JSON.stringify(report,null,2)+'\n');
console.log(JSON.stringify(report,null,2));
process.exitCode=report.passed?0:1;
