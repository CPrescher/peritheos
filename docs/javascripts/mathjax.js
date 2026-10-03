window.MathJax = {
  tex: {
    // Arithmatex wraps Markdown math; notebook exports retain dollar delimiters.
    inlineMath: [["\\(", "\\)"], ["$", "$"]],
    displayMath: [["\\[", "\\]"], ["$$", "$$"]],
    processEscapes: true,
    processEnvironments: true,
  },
  options: {
    // Allow traversal into notebook paragraphs, while excluding code and copies.
    ignoreHtmlClass: "highlight|jp-CodeCell|clipboard-copy-txt",
    processHtmlClass: "arithmatex|jp-RenderedMarkdown",
  },
};
