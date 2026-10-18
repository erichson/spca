#!/usr/bin/env node
// Optional development check; requires playwright-core and a Chromium browser.
// Usage: node check-github-rendering.cjs COMMIT README_SNAPSHOT OUTPUT_DIR [OWNER/REPO]
// OUTPUT_DIR must be absent or empty. This checks the public, immutable GitHub
// page and verifies that README_SNAPSHOT matches its raw source byte for byte.
// It uses a new browser context without an account, cookies or credentials.
const fs = require('fs');
const path = require('path');
const crypto = require('crypto');

async function launchBrowser() {
  let chromium;
  try {
    ({ chromium } = require('playwright-core'));
  } catch (_) {
    throw new Error('Install the optional playwright-core development tool, or set NODE_PATH to the directory containing it.');
  }
  return chromium.launch({
    executablePath: process.env.SPARSEPCA_CHROMIUM_EXECUTABLE || undefined,
    headless: true,
    args: ['--disable-dev-shm-usage']
  });
}

// This function is passed directly to Playwright and runs inside the page.
function collectDOM() {
  const article = document.querySelector('article.markdown-body');
  if (!article) return { missingArticle: true };
  // Word boundaries distinguish a renderer's ParseError/parse error from
  // ordinary prose such as "sparse error term S".
  const errors = /(?:following macros are not allowed|macro[^\n]{0,60}not allowed|unknown (?:command|macro)|undefined control sequence|\bparse\s*error\b|error (?:rendering|parsing) (?:math|latex)|(?:failed|unable) to render)/i;
  const errorSelector = '.flash-error, .katex-error, .MathJax_Error, merror, [data-mjx-error], [role="alert"]';
  const box = element => {
    const bounds = element.getBoundingClientRect();
    return { width: bounds.width, height: bounds.height };
  };
  const math = [...article.querySelectorAll('math-renderer')].map((element, index) => {
    const root = element.shadowRoot || element;
    const mathml = [...root.querySelectorAll('math')].map(m => ({
      text: m.textContent,
      tokens: [...m.querySelectorAll('mi,mn,mo,mtext')].map(t => t.textContent),
      ...box(m),
      display: getComputedStyle(m).display,
      visibility: getComputedStyle(m).visibility
    }));
    const svgs = [...root.querySelectorAll('svg')].map(svg => {
      const shapes = [...svg.querySelectorAll('path,rect,line,polyline,polygon,circle,ellipse,text,use')]
        .filter(n => !n.closest('defs'));
      const brokenUses = [...svg.querySelectorAll('use')]
        .map(u => u.getAttribute('href') || u.getAttribute('xlink:href'))
        .filter(href => !href || (href.startsWith('#') &&
          !root.querySelector('[id="' + CSS.escape(href.slice(1)) + '"]') &&
          !document.getElementById(href.slice(1))));
      let paintedBounds = null;
      try {
        const bounds = svg.getBBox();
        paintedBounds = { width: bounds.width, height: bounds.height };
      } catch (_) { /* A missing bounding box is reported by assess(). */ }
      return {
        ...box(svg), viewBox: svg.getAttribute('viewBox'), shapes: shapes.length,
        brokenUses, paintedBounds, display: getComputedStyle(svg).display,
        visibility: getComputedStyle(svg).visibility
      };
    });
    return {
      index, text: root.textContent, markup: root.innerHTML, ...box(element), mathml, svgs,
      errorElements: [...root.querySelectorAll(errorSelector)].map(n => ({
        tag: n.tagName, class: String(n.className), text: n.textContent
      })),
      errorText: errors.test(root.textContent),
      rawTexFallback: !!root.querySelector('pre,code')
    };
  });
  return {
    title: document.title, url: location.href, article: box(article), math,
    visibleErrorText: errors.test(article.innerText),
    errorElements: [...article.querySelectorAll(errorSelector)].map(n => ({
      tag: n.tagName, class: String(n.className), text: n.textContent
    })),
    images: [...article.querySelectorAll('img')].map(i => ({
      src: i.currentSrc, alt: i.alt, complete: i.complete,
      width: i.naturalWidth, height: i.naturalHeight
    })),
    unresolvedAnchors: [...article.querySelectorAll('a[href^="#"]')]
      .map(a => a.getAttribute('href')).filter(href => {
        if (href === '#') return false;
        const id = decodeURIComponent(href.slice(1));
        return !document.getElementById(id) && !document.getElementById('user-content-' + id);
      }),
    bodyText: article.innerText
  };
}

function formulas(markdown) {
  return [...markdown.matchAll(/^```math\s*\n([\s\S]*?)^```\s*$/gm)].map(m => m[1].trim());
}

function assess(data, source) {
  const problems = [];
  const expected = formulas(source);
  if (data.missingArticle) return ['No rendered README article found.'];
  if (data.math.length !== expected.length) {
    problems.push(`Expected ${expected.length} math blocks; found ${data.math.length} wrappers.`);
  }
  if (data.visibleErrorText) problems.push('Rendered article contains a math error or macro-rejection message.');
  for (const m of data.math) {
    if (m.errorText || m.errorElements.length) problems.push(`Equation ${m.index + 1}: error text or error element.`);
    if (m.rawTexFallback) problems.push(`Equation ${m.index + 1}: unrendered TeX fallback remains.`);
    const hasMathML = m.mathml.some(x => x.text.trim() && x.width > 0 && x.height > 0 &&
      x.display !== 'none' && x.visibility !== 'hidden');
    const hasSVG = m.svgs.some(s => s.shapes > 0 && s.width > 0 && s.height > 0 &&
      s.paintedBounds?.width > 0 && s.paintedBounds?.height > 0 && !s.brokenUses.length &&
      s.display !== 'none' && s.visibility !== 'hidden');
    if (!hasMathML && !hasSVG) problems.push(`Equation ${m.index + 1}: no nonempty visible MathML or painted SVG output.`);
    for (const svg of m.svgs) {
      if (!svg.shapes || svg.brokenUses.length || !(svg.paintedBounds?.width > 0) || !(svg.paintedBounds?.height > 0)) {
        problems.push(`Equation ${m.index + 1}: empty or broken SVG output.`);
      }
    }
    // The equation involved in the reported bug needs a semantic check as well
    // as a rendered box. NFKC normalizes MathML's styled Unicode characters.
    if (/soft/.test(expected[m.index] || '') && /sign/.test(expected[m.index] || '')) {
      const text = m.mathml.map(x => x.tokens.join('')).join('').normalize('NFKC')
        .replace(/[\s\u2061\u2062\u2063]/g, '');
      const tokens = ['soft', 'sign', 'max', 'u', 'τ', '0', '=', '|', '−'];
      const missing = tokens.filter(token => !text.includes(token));
      if (missing.length) problems.push(`Equation ${m.index + 1}: rendered soft-threshold formula lacks ${missing.join(', ')}.`);
      m.semanticCheck = { expectedTokens: tokens, normalizedMathML: text, missing };
    }
  }
  for (const image of data.images) {
    if (!image.complete || !image.width || !image.height) problems.push(`Image did not render: ${image.src}`);
  }
  if (data.unresolvedAnchors.length) problems.push('Unresolved anchors: ' + data.unresolvedAnchors.join(', '));
  return problems;
}

async function inspect(commit, readmePath, outputDir, repository = 'erichson/spca') {
  if (!/^[a-f0-9]{40}$/.test(commit)) throw new Error('Pass the full immutable source commit.');
  if (!/^[A-Za-z0-9_.-]+\/[A-Za-z0-9_.-]+$/.test(repository)) throw new Error('Repository must have the form OWNER/REPO.');
  const sourceBytes = fs.readFileSync(readmePath);
  const source = sourceBytes.toString('utf8');
  const out = path.resolve(outputDir);
  if (fs.existsSync(out) && fs.readdirSync(out).length) throw new Error('OUTPUT_DIR must be absent or empty to preserve earlier evidence.');
  fs.mkdirSync(out, { recursive: true });
  const ignoreHTTPSErrors = process.env.SPARSEPCA_PUBLIC_BROWSER_IGNORE_HTTPS_ERRORS === '1';
  const browser = await launchBrowser();
  try {
    const context = await browser.newContext({
      viewport: { width: 1440, height: 1000 }, deviceScaleFactor: 1,
      ignoreHTTPSErrors
    });
    const rawUrl = `https://raw.githubusercontent.com/${repository}/${commit}/README.md`;
    // Navigate with the browser itself so source verification uses the same
    // proxy and certificate configuration as the rendered-page inspection.
    const sourcePage = await context.newPage();
    const rawResponse = await sourcePage.goto(rawUrl, { waitUntil: 'domcontentloaded', timeout: 60000 });
    if (rawResponse?.status() !== 200) throw new Error(`README source returned HTTP ${rawResponse?.status()}: ${rawUrl}`);
    const rawBytes = await rawResponse.body();
    await sourcePage.close();
    const sha256 = bytes => crypto.createHash('sha256').update(bytes).digest('hex');
    const sourceVerification = {
      rawUrl, httpStatus: rawResponse.status(), bytes: rawBytes.length,
      localSHA256: sha256(sourceBytes), remoteSHA256: sha256(rawBytes),
      matched: sourceBytes.equals(rawBytes)
    };
    fs.writeFileSync(path.join(out, 'source-verification.json'), JSON.stringify(sourceVerification, null, 2) + '\n');
    if (!sourceVerification.matched) throw new Error('README_SNAPSHOT differs from the README at the requested GitHub commit.');
    const page = await context.newPage();
    const failedRequests = [], pageErrors = [];
    page.on('requestfailed', r => failedRequests.push({ url: r.url(), error: r.failure()?.errorText }));
    page.on('pageerror', e => pageErrors.push(e.message));
    const url = `https://github.com/${repository}/blob/${commit}/README.md`;
    const response = await page.goto(url, { waitUntil: 'domcontentloaded', timeout: 60000 });
    await page.locator('article.markdown-body').waitFor({ timeout: 30000 });
    await page.locator('article.markdown-body details').evaluateAll(nodes => nodes.forEach(n => { n.open = true; }));
    await page.waitForLoadState('networkidle', { timeout: 30000 }).catch(() => {});
    await page.evaluate(() => document.fonts.ready);
    // Wait for output or an explicit failure, then assess what actually exists.
    await page.waitForFunction(() => [...document.querySelectorAll('article.markdown-body math-renderer')].every(m => {
      const root = m.shadowRoot || m;
      return root.querySelector('math,svg,.flash-error,.katex-error,merror,pre');
    }), {}, { timeout: 30000 }).catch(() => {});
    const data = await page.evaluate(collectDOM);
    const problems = assess(data, source);
    if (response?.status() !== 200) problems.push(`HTTP status ${response?.status()}.`);
    if (failedRequests.length) problems.push(`${failedRequests.length} failed network request(s).`);
    if (pageErrors.length) problems.push(`${pageErrors.length} browser page error(s).`);
    fs.writeFileSync(path.join(out, 'github-readme.html'), await page.content());
    await page.screenshot({ path: path.join(out, 'github-desktop-full.png'), fullPage: true });
    const softIndex = formulas(source).findIndex(f => /soft/.test(f));
    if (softIndex >= 0 && data.math[softIndex]) {
      const equation = page.locator('article.markdown-body math-renderer').nth(softIndex);
      await equation.scrollIntoViewIfNeeded();
      await page.screenshot({ path: path.join(out, 'soft-threshold-context.png') });
      await equation.screenshot({ path: path.join(out, 'soft-threshold-equation.png') });
    }
    await page.setViewportSize({ width: 390, height: 844 });
    await page.evaluate(() => window.scrollTo(0, 0));
    await page.screenshot({ path: path.join(out, 'github-mobile-top.png') });
    const mobile = await page.evaluate(() => ({
      width: document.documentElement.clientWidth,
      scrollWidth: document.documentElement.scrollWidth,
      articleWidth: document.querySelector('article.markdown-body').clientWidth,
      articleScrollWidth: document.querySelector('article.markdown-body').scrollWidth
    }));
    if (mobile.scrollWidth > mobile.width) problems.push(`Mobile page overflow: ${mobile.scrollWidth} > ${mobile.width}.`);
    const result = {
      repository, requestedUrl: url, sourceCommit: commit, sourceVerification,
      httpStatus: response?.status(), browserVersion: browser.version(),
      publicUnauthenticatedContext: true, ignoreHTTPSErrors,
      certificateNote: ignoreHTTPSErrors ?
        'Certificate verification was explicitly disabled for this isolated public-page check; no signed-in session or credentials were used.' :
        'Certificate verification enabled.',
      inspectedAtUTC: new Date().toISOString(), expectedMath: formulas(source),
      failedRequests, pageErrors, mobile, ...data, problems, passed: !problems.length
    };
    fs.writeFileSync(path.join(out, 'github-inspection.json'), JSON.stringify(result, null, 2) + '\n');
    console.log(JSON.stringify({
      sourceCommit: commit, passed: result.passed, sourceMatched: sourceVerification.matched,
      mathWrappers: data.math.length, renderedMathML: data.math.filter(m => m.mathml.length).length,
      images: data.images.length, problems
    }, null, 2));
    return result.passed;
  } finally {
    await browser.close();
  }
}

module.exports = { launchBrowser, collectDOM, formulas, assess, inspect };
if (require.main === module) {
  const [commit, readme, out, repository] = process.argv.slice(2);
  if (!commit || !readme || !out) {
    console.error('Usage: node check-github-rendering.cjs COMMIT README_SNAPSHOT OUTPUT_DIR [OWNER/REPO]');
    process.exitCode = 2;
  } else {
    inspect(commit, readme, out, repository).then(ok => { if (!ok) process.exitCode = 1; })
      .catch(error => { console.error(error.message); process.exitCode = 2; });
  }
}
