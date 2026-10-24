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
const FormData = class { constructor(form) { this.form=form; } entries() { return this.form.data; } };
const context = vm.createContext({document, console, FormData});
const prefix = source.slice(0,source.indexOf('function renderPolytropicCharts('));
assert.ok(prefix.includes('function renderLineChart('));
vm.runInContext(prefix, context);
if (source.includes('function renderGeometryExportScope(')) {
  vm.runInContext(source.slice(source.indexOf('function renderGeometryExportScope('),
    source.indexOf('function renderGeometryComparison(')),context);
}
const renderSource = source.slice(source.indexOf('function renderUncertainty('),source.indexOf('function analysisWorkspace('));
assert.ok(renderSource.startsWith('function renderUncertainty('));
vm.runInContext('let latestUncertaintyResult=null;\n'+renderSource,context);
const cycleStart=source.includes('function updateCycleControls(')
  ? source.indexOf('function updateCycleControls(') : source.indexOf('function renderEngineCycle(');
vm.runInContext('let latestCycleResult=null;\n'+source.slice(cycleStart,source.indexOf('function renderUncertainty(')),context);
vm.runInContext(source.slice(source.indexOf('function renderLossBudget('),source.indexOf('function renderTwoPhaseBand(')),context);
vm.runInContext(source.slice(source.indexOf('function renderTwoPhaseBand('),source.indexOf('function renderNozzleThermal(')),context);
vm.runInContext(source.slice(source.indexOf('function renderNozzleThermal('),source.indexOf('let latestStagingResult')),context);
vm.runInContext('let latestStagingResult=null;let latestBurnResult=null;let latestBurnInput=null;\n'+
  source.slice(source.indexOf('function renderStaging('),source.indexOf('function updateCycleControls(')),context);
vm.runInContext(source.slice(source.indexOf('function analysisWorkspace('),source.indexOf('function downloadText(')),context);
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
check('loss covariance bands retain proxy qualification, basis and literal warning text',()=>{
  context.renderLossBudget({operating_basis:'fixed-pressure-area',mass_flow_ratio:.8,coverage_qualified:false,
    delivered_thrust_band_n:[900,1100],terms:[{thrust_loss_n:0},{thrust_loss_n:-10}],
    warnings:[{message:'<img onerror=bad>'}]});
  const note=document.getElementById('loss-scope-note');
  assert.match(note.textContent,/Pc.*At/);assert.match(note.textContent,/kapsama.*değildir/);
  assert.match(note.textContent,/<img onerror=bad>/);assert.equal(note.children.length,0);
});
check('legacy loss basis visibly distinguishes its multiplier from physical discharge',()=>{
  context.renderLossBudget({operating_basis:'fixed-mass-flow-lumped',mass_flow_ratio:1,
    delivered_thrust_band_n:[900,1100],terms:[{thrust_loss_n:1}],warnings:[]});
  assert.match(document.getElementById('loss-scope-note').textContent,/debi katsayısı değildir/);
});
check('loss form sends the correlation matrix as an array and refuses malformed JSON',()=>{
  const form={id:'loss-budget-form',data:[['operating_basis','fixed-pressure-area'],
    ['correlation_matrix','[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]']]};
  const data=context.formPayload(form);
  assert.ok(Array.isArray(data.correlation_matrix));assert.equal(data.correlation_matrix[2][2],1);
  form.data[1][1]='[[invalid]]';assert.throws(()=>context.formPayload(form));
});
check('two-phase energy means specific energy and keeps pressure separate from momentum',()=>{
  context.renderTwoPhaseBand({coupling_definition:'specific-kinetic-energy-ratio',reference_pressure_thrust_n:100,
    thrust_range_n:[975,1000],specific_impulse_range_s:[400,410],contains_nonpositive_thrust:false,warnings:[]});
  const note=document.getElementById('two-phase-scope-note');
  assert.match(note.textContent,/birim kütle/);assert.match(note.textContent,/100/);
  assert.match(note.textContent,/kapsama.*değildir/);
});
check('nonpositive two-phase endpoints and zero-pressure assumption stay visible as literal text',()=>{
  context.renderTwoPhaseBand({coupling_definition:'exit-velocity-ratio',reference_pressure_thrust_n:0,
    thrust_range_n:[-1,1],specific_impulse_range_s:[-450,450],contains_nonpositive_thrust:true,
    warnings:[{message:'<img onerror=bad>'}]});
  const note=document.getElementById('two-phase-scope-note');
  assert.match(note.textContent,/negatif/);assert.match(note.textContent,/varsayım/);
  assert.match(note.textContent,/<img onerror=bad>/);assert.equal(note.children.length,0);
});
check('staging distinguishes absent and non-nested grid comparisons from convergence evidence',()=>{
  context.renderStaging({stages:[{delta_v_m_s:100}],refinement_performed:false});
  assert.match(document.getElementById('staging-scope-note').textContent,/mevcut değil/);
  context.renderStaging({stages:[{delta_v_m_s:100}],refinement_performed:true,refinement_grid_nested:false});
  assert.match(document.getElementById('staging-scope-note').textContent,/iç içe değil/);
  assert.match(document.getElementById('staging-scope-note').textContent,/hata sınırı değildir/);
});
check('staging workspace retains solved inputs and full artifact after the form changes',()=>{
  vm.runInContext('latestCycleResult=null;latestUncertaintyResult=null;',context);
  document.getElementById('staging-form').data=[['target_delta_v_m_s','999']];
  const data={inputs:{target_delta_v_m_s:100,resolution:40},stages:[{delta_v_m_s:100}],
    initial_mass_kg:2,payload_fraction:.5,achieved_delta_v_m_s:100,active_constraints:[],
    refinement_performed:true,refinement_grid_nested:true,continuous_optimum_qualified:false};
  context.renderStaging(data);
  const saved=context.analysisWorkspace().cases[0];
  assert.equal(saved.inputs.target_delta_v_m_s,100);
  assert.equal(saved.results.artifact,data);
  assert.equal(saved.model_versions['staging-optimizer'],'2.0');
  assert.ok(saved.warnings.some(w=>w.includes('unqualified')));
});
check('thermal signed flux and unqualified ledger closure remain explicit with no extra warnings',()=>{
  context.renderNozzleThermal({stations:[{x_m:0,heat_flux_w_m2:1,adiabatic_wall_temperature_k:4},
    {x_m:1,heat_flux_w_m2:-1,adiabatic_wall_temperature_k:2}],contains_negative_heat_flux:true,warnings:[]});
  const note=document.getElementById('thermal-warning');
  assert.match(note.textContent,/negatif/);assert.match(note.textContent,/özdeşliğidir/);
  assert.match(note.textContent,/ölçülmüş enerji kapanışı değildir/);
  assert.match(note.textContent,/verilen duvar sıcaklığını/);
});
check('thermal warnings are literal text and never imply validated coupled wall cooling',()=>{
  context.renderNozzleThermal({stations:[{x_m:0,heat_flux_w_m2:1,adiabatic_wall_temperature_k:4},
    {x_m:1,heat_flux_w_m2:2,adiabatic_wall_temperature_k:3}],warnings:[{message:'<img onerror=bad>'}],
    cooling:{warnings:[],outlet_temperature_k:200,pressure_drop_pa:1,energy_balance_residual_w:0}});
  const note=document.getElementById('thermal-warning');
  assert.match(note.textContent,/<img onerror=bad>/);assert.equal(note.children.length,0);
  assert.match(note.textContent,/eşlenik.*yapılmaz/);
  assert.match(note.textContent,/yeterliliği doğrulanmadı/);
});
check('coolant exact limit flag stays critical when displayed outlet equals the supplied limit',()=>{
  context.renderNozzleThermal({stations:[{x_m:0,heat_flux_w_m2:1,adiabatic_wall_temperature_k:4},
    {x_m:1,heat_flux_w_m2:2,adiabatic_wall_temperature_k:3}],warnings:[],cooling:{
      outlet_temperature_k:1.0000000000000002,temperature_limit_exceeded:true,warnings:[],
      energy_balance_residual_w:0,displayed_temperature_energy_residual_w:1}});
  assert.match(document.getElementById('coolant-scope-note').textContent,/AŞILDI/);
  assert.match(document.getElementById('coolant-scope-note').textContent,/aynı yuvarlanabilir/);
  assert.equal(document.getElementById('coolant-energy-residual').textContent,'1');
});
check('coolant missing decisions are unknown and transition or shape assumptions remain explicit',()=>{
  const data={stations:[{x_m:0,heat_flux_w_m2:1,adiabatic_wall_temperature_k:4},
    {x_m:1,heat_flux_w_m2:2,adiabatic_wall_temperature_k:3}],warnings:[],cooling:{warnings:[]}};
  context.renderNozzleThermal(data);
  assert.match(document.getElementById('coolant-scope-note').textContent,/mevcut değil/);
  data.cooling.flow_regime='transition';context.renderNozzleThermal(data);
  assert.match(document.getElementById('coolant-scope-note').textContent,/doğrulanmış geçiş bağıntısı değildir/);
  data.cooling.flow_regime='laminar';data.cooling.inputs={laminar_darcy_poiseuille_number:64};
  context.renderNozzleThermal(data);
  assert.match(document.getElementById('coolant-scope-note').textContent,/dairesel eşdeğer varsayımı/);
});
check('geometry exports identify fluid domain, absent wall and STL metre import',()=>{
  context.renderGeometryExportScope({coordinate_unit:'m',manufacturing_qualified:false,physical_model_qualified:false});
  const note=document.getElementById('geometry-export-scope-note').textContent;
  assert.match(note,/kapalı akış hacmi/);assert.match(note,/Duvar kalınlığı veya üretim çizimi içermez/);
  assert.match(note,/STL içe aktarırken metre birimini seçin/);assert.match(note,/SVG.*ayrı ölçeklenir/);
  assert.match(note,/fiziksel doğrulaması açık/);
});
check('geometry missing units stay unknown even if untrusted qualification flags are true',()=>{
  context.renderGeometryExportScope({manufacturing_qualified:true,physical_model_qualified:true});
  const note=document.getElementById('geometry-export-scope-note').textContent;
  assert.match(note,/birimini doğrulayın/);assert.match(note,/fiziksel doğrulaması açık/);
  assert.doesNotMatch(note,/metre cinsindedir/);
});
const report={scope:'Actual chart/render/control function executions with instrumented DOM; no browser layout or physics claim',
  source_sha256:crypto.createHash('sha256').update(source).digest('hex'),node:process.version,rows,
  passed:rows.every(row=>row.status==='passed')};
const target=process.argv[2];
if (target) fs.writeFileSync(target,JSON.stringify(report,null,2)+'\n');
console.log(JSON.stringify(report,null,2));
process.exitCode=report.passed?0:1;
