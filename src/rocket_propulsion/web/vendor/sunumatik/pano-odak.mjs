/* pano-odak.mjs — "odak düzeni"nin davranışı (bkz. core/pano-odak.css).

   odakKur({ kok }) :
     · .po-sekmeler içindeki role=tab düğmelerini aria-controls ile panellere bağlar
       (tıklama, ←/→/Home/End; seçili sekme tabindex 0, diğerleri −1).
     · Sekme seçilince açık "Ayrıntı" kapanır (grafik alanını geri verir).
     · #caption (yöntem metni) varsa [data-po-yontem] içine taşınır: id ve içerik
       aynı kalır, sayfa onu eskisi gibi yazar.
     · "Daha fazla ayar" (details.po-ayar): Esc ve dışarı tıklama kapatır.
   Sayfa grafikleri panel görünür olunca boyutlanır: sayfalar ResizeObserver ile
   .po-panel öğelerini de izler (gizli panel 0 px, açılınca gerçek boyut). */

export function odakKur({ kok = document } = {}) {
  for (const liste of kok.querySelectorAll('.po-sekmeler')) {
    const sekmeler = [...liste.querySelectorAll('[role="tab"]')];
    const sec = (hedef, odakla = false) => {
      for (const s of sekmeler) {
        const acik = s === hedef;
        s.setAttribute('aria-selected', String(acik));
        s.tabIndex = acik ? 0 : -1;
        const panel = kok.getElementById ? kok.getElementById(s.getAttribute('aria-controls')) : document.getElementById(s.getAttribute('aria-controls'));
        if (panel) panel.hidden = !acik;
      }
      const ayrinti = liste.closest('.po-yan')?.querySelector('.po-ayrinti');
      if (ayrinti) ayrinti.open = false;
      if (odakla) hedef.focus();
    };
    const ilk = sekmeler.find(s => s.getAttribute('aria-selected') === 'true') ?? sekmeler[0];
    if (ilk) sec(ilk);
    for (const s of sekmeler) s.addEventListener('click', () => sec(s));
    liste.addEventListener('keydown', e => {
      const i = sekmeler.indexOf(document.activeElement); if (i < 0) return;
      const n = sekmeler.length;
      const j = e.key === 'ArrowRight' ? (i + 1) % n : e.key === 'ArrowLeft' ? (i - 1 + n) % n : e.key === 'Home' ? 0 : e.key === 'End' ? n - 1 : -1;
      if (j < 0) return;
      e.preventDefault(); e.stopPropagation(); sec(sekmeler[j], true);
    });
  }
  const yontem = kok.querySelector('[data-po-yontem]'), caption = kok.querySelector('#caption');
  if (yontem && caption && !yontem.contains(caption)) yontem.appendChild(caption);
  for (const ayar of kok.querySelectorAll('details.po-ayar')) {
    ayar.addEventListener('keydown', e => { if (e.key === 'Escape' && ayar.open) { e.stopPropagation(); ayar.open = false; ayar.querySelector('summary')?.focus(); } });
    document.addEventListener('pointerdown', e => { if (ayar.open && !ayar.contains(e.target)) ayar.open = false; });
  }
}
