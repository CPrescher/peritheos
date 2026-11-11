/* Preserve bookmarks to sections of the former single-page equation reference.
 * The index also retains a visible family link beside every legacy anchor.
 */
(() => {
  if (!window.location.hash) return;
  let anchor;
  try {
    anchor = decodeURIComponent(window.location.hash.slice(1));
  } catch {
    return;
  }
  const target = document.getElementById(anchor)?.dataset.equationTarget;
  if (target) window.location.replace(new URL(target, window.location.href));
})();
