"use strict";

const DATA = window.EC_TEMPLATE_DATA;
const TEMPLATE_LIMIT_STEP = 60;
const PREFIX_BROWSE_DEPTH = 3;

const state = {
  prefix: "",
  selectedEc: null,
  query: "",
  templateLimit: TEMPLATE_LIMIT_STEP,
};

const ecEntries = DATA.ecEntries;
const ecByCode = new Map(ecEntries.map((entry) => [entry.ec, entry]));
const templatesByEc = new Map();
for (const template of DATA.templates) {
  if (!templatesByEc.has(template.ec)) {
    templatesByEc.set(template.ec, []);
  }
  templatesByEc.get(template.ec).push(template);
}

const els = {};

document.addEventListener("DOMContentLoaded", init);

function init() {
  for (const id of [
    "datasetMeta",
    "rdkitStatus",
    "searchInput",
    "clearSearch",
    "resetBrowse",
    "breadcrumbs",
    "browseList",
    "panelHead",
    "resultSummary",
    "ecResults",
    "templateList",
    "loadMore",
    "backToTop",
  ]) {
    els[id] = document.getElementById(id);
  }

  els.datasetMeta.textContent = `Supported by ARIES DB: ${DATA.metadata.uniqueEcCount.toLocaleString()} EC buckets, ${DATA.metadata.uniqueTemplateCount.toLocaleString()} templates, ${DATA.metadata.sourceCsvRows.toLocaleString()} source rows`;

  els.searchInput.addEventListener("input", () => {
    state.query = els.searchInput.value.trim();
    state.templateLimit = TEMPLATE_LIMIT_STEP;
    render();
  });
  els.clearSearch.addEventListener("click", () => {
    els.searchInput.value = "";
    state.query = "";
    state.templateLimit = TEMPLATE_LIMIT_STEP;
    render();
  });
  els.resetBrowse.addEventListener("click", () => {
    state.prefix = "";
    state.selectedEc = null;
    state.templateLimit = TEMPLATE_LIMIT_STEP;
    render();
  });
  els.loadMore.addEventListener("click", () => {
    state.templateLimit += TEMPLATE_LIMIT_STEP;
    renderMainPanel();
  });
  els.backToTop.addEventListener("click", () => {
    window.scrollTo({ top: 0, behavior: "smooth" });
  });
  window.addEventListener("scroll", updateBackToTop);

  const hashEc = new URLSearchParams(window.location.hash.slice(1)).get("ec");
  if (hashEc && ecByCode.has(hashEc)) {
    selectEc(hashEc, { clearSearch: true });
  }

  window.rdkitReady
    .then((RDKit) => {
      els.rdkitStatus.textContent = `RDKit.js ${RDKit.version()}`;
      els.rdkitStatus.classList.add("ready");
      renderVisibleDrawings();
    })
    .catch(() => {
      els.rdkitStatus.textContent = "RDKit.js unavailable";
      els.rdkitStatus.classList.add("error");
    });

  render();
  updateBackToTop();
}

function render() {
  renderBrowse();
  renderMainPanel();
}

function renderBrowse() {
  renderBreadcrumbs();
  const nodes = getBrowseNodes(state.prefix);
  els.browseList.textContent = "";

  for (const node of nodes) {
    const button = document.createElement("button");
    button.type = "button";
    button.className = "browse-item";
    if (node.code === state.selectedEc || node.code === state.prefix) {
      button.classList.add("active");
    }
    button.innerHTML = `
      <span>
        <span class="browse-code">${escapeHtml(node.code)}</span>
        <span class="browse-name">${escapeHtml(node.name || "No name")}</span>
      </span>
      <span class="count-stack">
        <span>${node.ecCount.toLocaleString()} EC</span><br>
        <span>${node.templateCount.toLocaleString()} tpl</span>
      </span>
    `;
    button.addEventListener("click", () => {
      if (node.leaf) {
        selectEc(node.code, { clearSearch: true });
      } else {
        state.prefix = node.code;
        state.selectedEc = null;
        state.templateLimit = TEMPLATE_LIMIT_STEP;
        render();
      }
    });
    els.browseList.appendChild(button);
  }
}

function renderBreadcrumbs() {
  els.breadcrumbs.textContent = "";
  addCrumb("Top", "");
  if (!state.prefix) {
    return;
  }
  const parts = state.prefix.split(".");
  let current = "";
  for (const part of parts) {
    current = current ? `${current}.${part}` : part;
    addCrumb(current, current);
  }
}

function addCrumb(label, prefix) {
  const button = document.createElement("button");
  button.type = "button";
  button.textContent = label;
  button.addEventListener("click", () => {
    state.prefix = prefix;
    state.selectedEc = null;
    state.templateLimit = TEMPLATE_LIMIT_STEP;
    render();
  });
  els.breadcrumbs.appendChild(button);
}

function getBrowseNodes(prefix) {
  const nodes = new Map();
  for (const entry of ecEntries) {
    const parts = splitEc(entry.ec);
    if (!parts.length) {
      if (!prefix) {
        nodes.set(entry.ec, leafNode(entry));
      }
      continue;
    }

    const prefixParts = prefix ? prefix.split(".") : [];
    if (prefix && parts.slice(0, prefixParts.length).join(".") !== prefix) {
      continue;
    }

    if (!prefix) {
      const code = parts[0];
      setPrefixNode(nodes, code);
      continue;
    }

    const next = parts[prefixParts.length];
    if (!next || prefixParts.length >= PREFIX_BROWSE_DEPTH || !isDigits(next)) {
      nodes.set(entry.ec, leafNode(entry));
    } else {
      setPrefixNode(nodes, `${prefix}.${next}`);
    }
  }
  return Array.from(nodes.values()).sort((a, b) => compareEc(a.code, b.code));
}

function splitEc(ec) {
  const parts = ec.split(".").filter(Boolean);
  return parts.length && isDigits(parts[0]) ? parts : [];
}

function setPrefixNode(nodes, code) {
  const prefix = DATA.prefixes[code];
  if (!prefix) {
    return;
  }
  nodes.set(code, {
    code,
    name: prefix.name,
    templateCount: prefix.templateCount,
    rowCount: prefix.rowCount,
    ecCount: prefix.ecCount,
    leaf: false,
  });
}

function leafNode(entry) {
  return {
    code: entry.ec,
    name: entry.name,
    templateCount: entry.templateCount,
    rowCount: entry.rowCount,
    ecCount: 1,
    leaf: true,
  };
}

function renderMainPanel() {
  const query = state.query.toLowerCase();
  els.ecResults.textContent = "";

  if (query) {
    renderSearch(query);
    return;
  }

  if (state.selectedEc) {
    const entry = ecByCode.get(state.selectedEc);
    const templates = templatesByEc.get(state.selectedEc) || [];
    renderEcPanel(entry, templates);
    renderTemplates(templates);
    return;
  }

  const prefix = state.prefix ? DATA.prefixes[state.prefix] : null;
  els.panelHead.innerHTML = `
    <div class="panel-kicker">Browse</div>
    <div class="panel-title-row">
      <div>
        <h3>${escapeHtml(state.prefix || "Top-level EC classes")}</h3>
        <div class="panel-name">${escapeHtml(prefix ? prefix.name : "Select an EC class or search the template set.")}</div>
      </div>
    </div>
  `;
  els.resultSummary.textContent = "";
  els.templateList.innerHTML = `<div class="empty">No EC number selected.</div>`;
  els.loadMore.parentElement.classList.remove("visible");
}

function renderSearch(query) {
  const ecMatches = ecEntries
    .filter((entry) => matchesText(query, [entry.ec, entry.name]))
    .sort((a, b) => compareEc(a.ec, b.ec));
  const templateMatches = DATA.templates.filter((template) =>
    matchesText(query, [
      template.ec,
      ecByCode.get(template.ec)?.name || "",
      template.templateId,
      template.siteSmarts,
      template.radius0Smarts,
    ]),
  );

  els.panelHead.innerHTML = `
    <div class="panel-kicker">Search</div>
    <div class="panel-title-row">
      <div>
        <h3>${escapeHtml(state.query)}</h3>
        <div class="panel-name">${ecMatches.length.toLocaleString()} EC matches, ${templateMatches.length.toLocaleString()} template matches</div>
      </div>
    </div>
  `;
  els.resultSummary.textContent = "";
  renderEcMatches(ecMatches.slice(0, 24));
  renderTemplates(templateMatches);
}

function renderEcMatches(entries) {
  els.ecResults.textContent = "";
  for (const entry of entries) {
    const button = document.createElement("button");
    button.type = "button";
    button.className = "ec-result";
    button.innerHTML = `
      <span class="ec-code">${escapeHtml(entry.ec)}</span>
      <span class="browse-name">${escapeHtml(entry.name)}</span>
      <span class="badge-row">
        <span class="badge">${entry.templateCount.toLocaleString()} templates</span>
        <span class="badge">${entry.rowCount.toLocaleString()} rows</span>
      </span>
    `;
    button.addEventListener("click", () => selectEc(entry.ec, { clearSearch: true }));
    els.ecResults.appendChild(button);
  }
}

function renderEcPanel(entry, templates) {
  const sourceLabel = nameSourceLabel(entry.nameSource);
  els.panelHead.innerHTML = `
    <div class="panel-kicker">EC number</div>
    <div class="panel-title-row">
      <div>
        <h3 class="ec-code">${escapeHtml(entry.ec)}</h3>
        <div class="panel-name">${escapeHtml(entry.name)}</div>
      </div>
      <button type="button" id="copyEc">Copy EC</button>
    </div>
    <div class="badge-row">
      <span class="badge">${templates.length.toLocaleString()} templates</span>
      <span class="badge">${entry.rowCount.toLocaleString()} source rows</span>
      <span class="badge ${entry.nameSource === "expasy-entry" ? "" : "warn"}">${escapeHtml(sourceLabel)}</span>
    </div>
  `;
  document.getElementById("copyEc").addEventListener("click", (event) => {
    copyText(entry.ec, event.currentTarget);
  });
  els.resultSummary.textContent = "";
}

function renderTemplates(templates) {
  const visible = templates.slice(0, state.templateLimit);
  els.templateList.textContent = "";

  if (!templates.length) {
    els.templateList.innerHTML = `<div class="empty">No matching templates.</div>`;
    els.loadMore.parentElement.classList.remove("visible");
    return;
  }

  els.resultSummary.textContent =
    visible.length === templates.length
      ? `${templates.length.toLocaleString()} templates shown`
      : `${visible.length.toLocaleString()} of ${templates.length.toLocaleString()} templates shown`;

  visible.forEach((template, index) => {
    els.templateList.appendChild(templateCard(template, index));
  });

  if (visible.length < templates.length) {
    els.loadMore.parentElement.classList.add("visible");
  } else {
    els.loadMore.parentElement.classList.remove("visible");
  }
  renderVisibleDrawings();
}

function templateCard(template, index) {
  const card = document.createElement("article");
  card.className = "template-card";
  card.dataset.templateIndex = String(index);
  const ecEntry = ecByCode.get(template.ec);
  const datasetText = Object.entries(template.sourceDatasets)
    .map(([name, count]) => `${name}: ${count}`)
    .join(", ");
  const showRadius = template.radius0Smarts !== template.siteSmarts;
  const examples = template.examples || [];
  const fullReactionCount = examples.filter((example) => example.fullReactionSmiles).length;

  card.innerHTML = `
    <div class="template-header">
      <div class="template-title-row">
        <div>
          <div class="template-id">${escapeHtml(template.templateId)}</div>
          <div class="browse-name"><span class="ec-code">${escapeHtml(template.ec)}</span> ${escapeHtml(ecEntry?.name || "")}</div>
        </div>
        <button type="button" class="copy-template" data-copy="${escapeAttribute(template.siteSmarts)}">Copy SMARTS</button>
      </div>
      <div class="badge-row">
        <span class="badge">${template.rowCount.toLocaleString()} rows</span>
        <span class="badge">${examples.length.toLocaleString()} examples</span>
        <span class="badge">${fullReactionCount.toLocaleString()} full reactions</span>
        <span class="badge">${escapeHtml(datasetText || "dataset")}</span>
        ${template.selectivityIssueCount ? `<span class="badge warn">${template.selectivityIssueCount.toLocaleString()} selectivity issue rows</span>` : ""}
      </div>
    </div>
    <div class="template-body">
      <div class="reaction-art" data-smarts="${escapeAttribute(template.siteSmarts)}">RDKit.js loading</div>
      ${smartsBlock("Reaction SMARTS", template.siteSmarts)}
      ${showRadius ? smartsBlock("Radius-0 SMARTS", template.radius0Smarts) : ""}
      <div class="detail-grid">
        <span><strong>Legal sites</strong> ${formatMaybe(template.legalSiteCount)}</span>
        <span><strong>Atoms</strong> ${formatMaybe(template.numAtoms)}</span>
        <span><strong>Positive atoms</strong> ${formatMaybe(template.positiveAtomCount)}</span>
        <span><strong>Legal atoms</strong> ${formatMaybe(template.legalAtomCount)}</span>
        <span><strong>Example IDs</strong> ${escapeHtml(template.exampleIds.join(", ") || "none")}</span>
      </div>
      ${examplesSection(template)}
    </div>
  `;

  attachCopyHandlers(card);
  const showMore = card.querySelector(".show-more-examples");
  if (showMore) {
    showMore.addEventListener("click", () => {
      const list = card.querySelector(".examples-list");
      const hiddenExamples = examples.slice(1);
      hiddenExamples.forEach((example, hiddenIndex) => {
        list.insertAdjacentHTML("beforeend", exampleCard(example, hiddenIndex + 2));
      });
      attachCopyHandlers(list);
      showMore.remove();
      renderVisibleDrawings();
    });
  }

  return card;
}

function attachCopyHandlers(root) {
  for (const button of root.querySelectorAll("[data-copy]")) {
    if (button.dataset.copyBound === "1") {
      continue;
    }
    button.dataset.copyBound = "1";
    button.addEventListener("click", (event) => {
      copyText(event.currentTarget.dataset.copy, event.currentTarget);
    });
  }
}

function examplesSection(template) {
  const examples = template.examples || [];
  if (!examples.length) {
    return `
      <section class="examples-block">
        <div class="examples-head">
          <span>Examples</span>
        </div>
        <div class="empty">No model-ready example metadata available for this template.</div>
      </section>
    `;
  }
  const remaining = examples.length - 1;
  return `
    <section class="examples-block">
      <div class="examples-head">
        <span>Examples</span>
        <span>${examples.length.toLocaleString()} source examples</span>
      </div>
      <div class="examples-list">
        ${exampleCard(examples[0], 1)}
      </div>
      ${
        remaining > 0
          ? `<button type="button" class="show-more-examples">Show ${remaining.toLocaleString()} more example${remaining === 1 ? "" : "s"}</button>`
          : ""
      }
    </section>
  `;
}

function exampleCard(example, index) {
  const reaction = example.fullReactionSmiles || "";
  const substrate = example.canonicalSubstrateSmiles || example.mappedSubstrateSmiles || "";
  const mappedSubstrate = example.mappedSubstrateSmiles || "";
  const hasDistinctMappedSubstrate = mappedSubstrate && mappedSubstrate !== substrate;
  const uniprotIds = Array.isArray(example.uniprotIds) ? example.uniprotIds.filter(Boolean) : [];
  return `
    <article class="example-card">
      <div class="example-topline">
        <strong>Example ${index}</strong>
        <span class="badge">${escapeHtml(example.dataset || "source")}</span>
        <span class="badge">${escapeHtml(example.id)}</span>
      </div>
      <div class="example-protein">
        <span class="example-label">UniProt ID</span>
        <span class="uniprot-links">${uniprotIdMarkup(uniprotIds)}</span>
      </div>
      <div class="example-grid">
        <div class="sequence-block">
          <div class="smarts-head">
            <span>Protein sequence</span>
            <button type="button" data-copy="${escapeAttribute(example.proteinSequence || "")}">Copy</button>
          </div>
          <pre><code>${escapeHtml(example.proteinSequence || "No sequence available")}</code></pre>
        </div>
        <div class="example-facts">
          <span>source pair: ${escapeHtml((example.sourcePairIds || []).join(", ") || "n/a")}</span>
          <span>reaction: ${escapeHtml((example.sourceReactionIds || []).join(", ") || "n/a")}</span>
          <span>enzyme: ${escapeHtml((example.sourceEnzymeIds || []).join(", ") || "n/a")}</span>
          <span>sequence length: ${formatMaybe((example.proteinSequence || "").length || null)}</span>
          <span>full reaction: ${reaction ? "available" : "not recovered"}</span>
        </div>
      </div>
      <div class="example-substrate-grid">
        ${compactValueBlock("Substrate SMILES", substrate)}
        ${hasDistinctMappedSubstrate ? compactValueBlock("Mapped substrate SMILES", mappedSubstrate) : ""}
      </div>
      ${
        reaction
          ? `
            <div class="reaction-art full-reaction-art" data-reaction="${escapeAttribute(reaction)}">RDKit.js loading</div>
            ${smartsBlock("Full mapped reaction SMILES", reaction)}
          `
          : `
            <div class="empty">Full mapped reaction was not available in the recovered source rows for this example.</div>
          `
      }
    </article>
  `;
}

function uniprotIdMarkup(ids) {
  if (!ids.length) {
    return `<span class="missing-value">Not available</span>`;
  }
  return ids
    .map((id) => {
      const text = String(id).trim();
      const url = `https://www.uniprot.org/uniprotkb/${encodeURIComponent(text)}/entry`;
      return `<a class="uniprot-link" href="${escapeAttribute(url)}" target="_blank" rel="noopener noreferrer">${escapeHtml(text)}</a>`;
    })
    .join("");
}

function compactValueBlock(label, value) {
  if (!value) {
    return "";
  }
  return `
    <div class="compact-value">
      <div class="smarts-head">
        <span>${escapeHtml(label)}</span>
        <button type="button" data-copy="${escapeAttribute(value)}">Copy</button>
      </div>
      <pre><code>${escapeHtml(value)}</code></pre>
    </div>
  `;
}

function smartsBlock(label, value) {
  return `
    <div class="smarts-block">
      <div class="smarts-head">
        <span>${escapeHtml(label)}</span>
        <button type="button" data-copy="${escapeAttribute(value)}">Copy</button>
      </div>
      <pre><code>${escapeHtml(value)}</code></pre>
    </div>
  `;
}

function renderVisibleDrawings() {
  if (!window.RDKit) {
    return;
  }
  for (const target of document.querySelectorAll(".reaction-art[data-smarts]")) {
    if (target.dataset.rendered === "1") {
      continue;
    }
    target.dataset.rendered = "1";
    renderReactionSmarts(target.dataset.smarts, target, "query");
  }
  for (const target of document.querySelectorAll(".reaction-art[data-reaction]")) {
    if (target.dataset.rendered === "1") {
      continue;
    }
    target.dataset.rendered = "1";
    renderReactionSmarts(target.dataset.reaction, target, "smiles");
  }
}

function renderReactionSmarts(smarts, target, mode = "query") {
  target.textContent = "";
  const parsed = parseReactionSmarts(smarts);
  if (!parsed) {
    target.innerHTML = `<div class="fallback-mol">${escapeHtml(smarts)}</div>`;
    return;
  }

  const row = document.createElement("div");
  row.className = "reaction-row";
  addMoleculeSide(row, parsed.reactants, mode);
  addOperator(row, "->");
  addMoleculeSide(row, parsed.products, mode);
  target.appendChild(row);
}

function addMoleculeSide(row, parts, mode) {
  parts.forEach((part, index) => {
    if (index > 0) {
      addOperator(row, "+");
    }
    row.appendChild(renderMolecule(part, mode));
  });
}

function addOperator(row, text) {
  const span = document.createElement("span");
  span.className = "operator";
  span.textContent = text;
  row.appendChild(span);
}

function renderMolecule(smarts, mode = "query") {
  const cell = document.createElement("div");
  cell.className = "mol-cell";
  let mol = null;
  try {
    mol =
      mode === "smiles"
        ? window.RDKit.get_mol(smarts) || window.RDKit.get_qmol(smarts)
        : window.RDKit.get_qmol(smarts) || window.RDKit.get_mol(smarts);
    if (mol) {
      cell.innerHTML = mol.get_svg(230, 160);
      mol.delete();
      return cell;
    }
  } catch (error) {
    if (mol && typeof mol.delete === "function") {
      mol.delete();
    }
  }
  const fallback = document.createElement("div");
  fallback.className = "fallback-mol";
  fallback.textContent = smarts;
  cell.appendChild(fallback);
  return cell;
}

function parseReactionSmarts(smarts) {
  let parts = smarts.split(">>");
  if (parts.length === 2) {
    return {
      reactants: splitSide(parts[0]),
      products: splitSide(parts[1]),
    };
  }
  parts = smarts.split(">");
  if (parts.length === 3) {
    return {
      reactants: splitSide(parts[0]),
      products: splitSide(parts[2]),
    };
  }
  return null;
}

function splitSide(side) {
  const text = stripOuterParens(side.trim());
  const parts = [];
  let start = 0;
  let parenDepth = 0;
  let bracketDepth = 0;
  for (let i = 0; i < text.length; i += 1) {
    const char = text[i];
    if (char === "[" && bracketDepth === 0) {
      bracketDepth += 1;
    } else if (char === "]" && bracketDepth > 0) {
      bracketDepth -= 1;
    } else if (!bracketDepth && char === "(") {
      parenDepth += 1;
    } else if (!bracketDepth && char === ")" && parenDepth > 0) {
      parenDepth -= 1;
    } else if (!bracketDepth && !parenDepth && char === ".") {
      parts.push(stripOuterParens(text.slice(start, i).trim()));
      start = i + 1;
    }
  }
  parts.push(stripOuterParens(text.slice(start).trim()));
  return parts.filter(Boolean);
}

function stripOuterParens(value) {
  let text = value;
  while (text.startsWith("(") && text.endsWith(")") && enclosesWholeString(text)) {
    text = text.slice(1, -1).trim();
  }
  return text;
}

function enclosesWholeString(text) {
  let depth = 0;
  let bracketDepth = 0;
  for (let i = 0; i < text.length; i += 1) {
    const char = text[i];
    if (char === "[") {
      bracketDepth += 1;
    } else if (char === "]" && bracketDepth > 0) {
      bracketDepth -= 1;
    } else if (!bracketDepth && char === "(") {
      depth += 1;
    } else if (!bracketDepth && char === ")") {
      depth -= 1;
      if (depth === 0 && i < text.length - 1) {
        return false;
      }
    }
  }
  return depth === 0;
}

function selectEc(ec, options = {}) {
  state.selectedEc = ec;
  state.prefix = parentPrefixForEc(ec);
  state.templateLimit = TEMPLATE_LIMIT_STEP;
  if (options.clearSearch) {
    state.query = "";
    els.searchInput.value = "";
  }
  window.history.replaceState(null, "", `#ec=${encodeURIComponent(ec)}`);
  render();
}

function parentPrefixForEc(ec) {
  const parts = splitEc(ec);
  if (parts.length >= 4) {
    return parts.slice(0, 3).join(".");
  }
  if (parts.length > 1) {
    return parts.slice(0, -1).join(".");
  }
  return "";
}

function matchesText(query, values) {
  return values.some((value) => String(value).toLowerCase().includes(query));
}

function compareEc(a, b) {
  const aParts = a.split(".");
  const bParts = b.split(".");
  const len = Math.max(aParts.length, bParts.length);
  for (let i = 0; i < len; i += 1) {
    const av = aParts[i] || "";
    const bv = bParts[i] || "";
    const aNum = isDigits(av);
    const bNum = isDigits(bv);
    if (aNum && bNum && Number(av) !== Number(bv)) {
      return Number(av) - Number(bv);
    }
    if (aNum !== bNum) {
      return aNum ? -1 : 1;
    }
    if (av !== bv) {
      return av.localeCompare(bv);
    }
  }
  return 0;
}

function isDigits(value) {
  return /^\d+$/.test(value);
}

function nameSourceLabel(source) {
  if (source === "expasy-entry") {
    return "Expasy entry";
  }
  if (source === "expasy-parent-class") {
    return "Expasy parent class";
  }
  if (source === "expasy-class") {
    return "Expasy class";
  }
  return "Dataset bucket";
}

function formatMaybe(value) {
  return value === null || value === undefined ? "n/a" : String(value);
}

async function copyText(text, button) {
  const original = button.textContent;
  try {
    if (navigator.clipboard?.writeText) {
      await navigator.clipboard.writeText(text);
    } else {
      fallbackCopy(text);
    }
    button.textContent = "Copied";
  } catch (error) {
    fallbackCopy(text);
    button.textContent = "Copied";
  }
  window.setTimeout(() => {
    button.textContent = original;
  }, 1000);
}

function fallbackCopy(text) {
  const textarea = document.createElement("textarea");
  textarea.value = text;
  textarea.setAttribute("readonly", "");
  textarea.style.position = "fixed";
  textarea.style.top = "-1000px";
  document.body.appendChild(textarea);
  textarea.select();
  document.execCommand("copy");
  textarea.remove();
}

function updateBackToTop() {
  if (!els.backToTop) {
    return;
  }
  els.backToTop.classList.toggle("visible", window.scrollY > 480);
}

function escapeHtml(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}

function escapeAttribute(value) {
  return escapeHtml(value).replaceAll("\n", "&#10;");
}
