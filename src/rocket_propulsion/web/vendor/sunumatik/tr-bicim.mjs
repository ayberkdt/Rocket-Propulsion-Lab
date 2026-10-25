/* tr-bicim.mjs — Türkçe sayfalarda GÖRÜNTÜLENEN sayıların ondalık ayırıcısı virgüldür.

   Sayfalar sayıları çoğunlukla `toFixed` ile yazıyordu (0.67 m/s, ε 0.85, L/D 0.3); aynı sayfada
   Intl tr-TR ile yazılanlar virgüllüydü — tek ekranda iki biçim. Binlerce çağrı noktasını tek tek
   değiştirmek hem eksik kalır hem de SVG koordinatı/CSS gibi NOKTA OLMASI GEREKEN dizeleri bozar.
   Bu modül yalnız insanın okuduğu katmana dokunur:
     · metin düğümleri (HTML ve SVG <text>; MutationObserver ile sonradan gelenler de),
     · tuval yazısı (fillText / strokeText; hizalama tutarlı kalsın diye measureText de).
   Öznitelikler (x="12.5", style), girdi değerleri, <code>/<pre>/<kbd>/<samp>/<textarea>/<script>/<style>
   ve [data-ham-sayi] taşıyan alt ağaçlar DOKUNULMAZ.

   Kural (trSayi): `\d+\.\d+` bir ondalıktır ve virgüle döner; ŞU İSTİSNALARLA:
     · binlik grup: 1–3 haneli tam kısım + nokta + TAM 3 hane (11.000, 924.647, 1.496.559) — Türkçe
       binlik ayırıcı noktadır (Intl tr-TR böyle yazar). "0.143" binlik olamaz → ondalık.
       Belirsizlik: toFixed(3) ile yazılmış 1.725 de bu biçimdedir — o çağrılar kaynakta tr biçimine
       çevrilir (scripts/validate-dil.mjs bunu denetler).
     · harf/nokta/virgül/eğik çizgi/# ile bitişik sayılar (v1.2, x.mjs, 1.2.3, 12/3.4, #1.5, 1.5em);
       üstel gösterimin mantisi ise virgül alır (1.2e-3 → 1,2e-3).
   Yalnız <html lang="tr"> sayfalarında etkindir; modül yüklenince bir kez kurulur. */

/* eksi işareti, cümle sonu noktası ve üstel gösterim (1.2e-3 → 1,2e-3) serbest */
const ONDALIK = /(^|[^\p{L}\p{N}_.,/#\\])(\d+)\.(\d+)(?!(?![eE][+\-−]?\d)[\p{L}\p{N}_]|\.\d)/gu;
const binlikMi = (tam, kesir) => kesir.length === 3 && tam.length <= 3 && tam[0] !== '0';

/** Bir dizedeki ondalık noktaları virgüle çevirir (binlik grupları korur). Saf. */
export function trSayi(s) {
  if (typeof s !== 'string' || s.indexOf('.') < 0) return s;
  return s.replace(ONDALIK, (m, on, tam, kesir) => (binlikMi(tam, kesir) ? m : `${on}${tam},${kesir}`));
}

const HARIC = 'code,pre,kbd,samp,textarea,script,style,[data-ham-sayi]';

function metinDuzelt(node) {
  const p = node.parentElement;
  if (!p || p.closest(HARIC)) return;
  const v = node.data, y = trSayi(v);
  if (y !== v) node.data = y;
}
function agacDuzelt(kok) {
  if (kok.nodeType === 3) { metinDuzelt(kok); return; }
  if (kok.nodeType !== 1 && kok.nodeType !== 9 && kok.nodeType !== 11) return;
  if (kok.nodeType === 1 && kok.closest(HARIC)) return;
  const tw = document.createTreeWalker(kok, NodeFilter.SHOW_TEXT);
  for (let n = tw.nextNode(); n; n = tw.nextNode()) metinDuzelt(n);
}

let kuruldu = false;
export function trBicimKur() {
  if (kuruldu || typeof document === 'undefined') return;
  kuruldu = true;
  /* tuval: metin çizen ve ölçen üç yöntem AYNI dönüşümü görür (ölçülen genişlik çizilenle eşleşsin) */
  for (const P of [globalThis.CanvasRenderingContext2D?.prototype, globalThis.OffscreenCanvasRenderingContext2D?.prototype]) {
    if (!P || P.__trBicim) continue;
    P.__trBicim = true;
    for (const ad of ['fillText', 'strokeText', 'measureText']) {
      const ozgun = P[ad];
      if (typeof ozgun !== 'function') continue;
      P[ad] = function (metin, ...geri) { return ozgun.call(this, typeof metin === 'string' ? trSayi(metin) : metin, ...geri); };
    }
  }
  const basla = () => {
    agacDuzelt(document.body || document.documentElement);
    new MutationObserver((kayitlar) => {
      for (const k of kayitlar) {
        if (k.type === 'characterData') metinDuzelt(k.target);
        else for (const n of k.addedNodes) agacDuzelt(n);
      }
    }).observe(document.documentElement, { subtree: true, childList: true, characterData: true });
  };
  if (document.body) basla(); else document.addEventListener('DOMContentLoaded', basla, { once: true });
}

if ((globalThis.document?.documentElement?.lang || '').toLowerCase().startsWith('tr')) trBicimKur();
