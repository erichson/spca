#!/usr/bin/env node
// Offline browser fixtures for GitHub's previously missed sanitizer error.
// Usage: node test-github-rendering.cjs [RESULTS_JSON]
const fs = require('fs');
const assert = require('assert');
const { launchBrowser, collectDOM, assess } = require('./check-github-rendering.cjs');

const source = '```math\n\\mathrm{soft}(u,\\tau)=\\mathrm{sign}(u)\\max(|u|-\\tau,0).\n```\n';
const correctFormula = '<math xmlns="http://www.w3.org/1998/Math/MathML" display="block">' +
  '<mi>soft</mi><mo>(</mo><mi>u</mi><mo>,</mo><mi>τ</mi><mo>)</mo><mo>=</mo>' +
  '<mi>sign</mi><mo>(</mo><mi>u</mi><mo>)</mo><mi>max</mi><mo>(</mo>' +
  '<mo>|</mo><mi>u</mi><mo>|</mo><mo>−</mo><mi>τ</mi><mo>,</mo><mn>0</mn><mo>)</mo></math>';
// GitHub produced this markup for README.md at commit
// 2996e81fa54e26ca4013d6c170175c0b09861de2. A math-renderer wrapper exists even
// though its content is this error and fallback, rather than an equation.
const rejectedFormula = '<!----><div class="flash flash-error"><!---->' +
  'The following macros are not allowed: operatorname<!----></div>\n' +
  '        <br>\n        <pre><!---->$$\\operatorname{soft}(u,\\tau)\n' +
  '= \\operatorname{sign}(u)\\max(|u|-\\tau,0).$$<!----></pre><!---->';

async function check() {
  const browser = await launchBrowser();
  const results = [];
  try {
    const page = await browser.newPage();
    // Fixtures are self-contained. Prevent accidental network access.
    await page.route('**/*', route => route.abort());
    const cases = [
      { name: 'Actual GitHub sanitizer failure', html: rejectedFormula, pass: false, expected: /error|unrendered|no nonempty/ },
      { name: 'Empty math wrapper', html: '', pass: false, expected: /no nonempty/ },
      { name: 'Empty SVG wrapper', html: '<svg width="100" height="30" viewBox="0 0 100 30"></svg>', pass: false, expected: /empty or broken SVG/ },
      { name: 'Broken SVG use reference', html: '<svg width="100" height="30" viewBox="0 0 100 30"><use href="#missing-glyph" /></svg>', pass: false, expected: /empty or broken SVG/ },
      { name: 'Visible formula with wrong semantics', html: '<math><mi>x</mi><mo>=</mo><mn>1</mn></math>', pass: false, expected: /formula lacks/ },
      { name: 'Correct rendered soft-threshold formula', html: correctFormula, pass: true },
      { name: 'Ordinary sparse error prose', prose: 'The sparse error term S absorbs large residuals.', html: correctFormula, pass: true },
      { name: 'ParseError message', prose: 'ParseError: invalid expression.', html: correctFormula, pass: false, expected: /math error or macro-rejection/ },
      { name: 'Spaced parse error message', prose: 'Parse error: invalid expression.', html: correctFormula, pass: false, expected: /math error or macro-rejection/ }
    ];
    for (const fixture of cases) {
      await page.setContent('<article class="markdown-body"><p>' + (fixture.prose || '') +
        '</p><math-renderer>' + fixture.html + '</math-renderer></article>');
      const data = await page.evaluate(collectDOM);
      const problems = assess(data, source);
      assert.strictEqual(!problems.length, fixture.pass, fixture.name + ': ' + problems.join('; '));
      if (fixture.expected) {
        assert(fixture.expected.test(problems.join('\n')), fixture.name + ': expected failure was not detected');
      }
      results.push({ name: fixture.name, expectedPass: fixture.pass, observedPass: !problems.length, problems });
    }
  } finally {
    await browser.close();
  }
  const report = { cases: results.length, passed: results.length, results };
  if (process.argv[2]) fs.writeFileSync(process.argv[2], JSON.stringify(report, null, 2) + '\n');
  console.log(JSON.stringify(report, null, 2));
}

check().catch(error => { console.error(error.message); process.exitCode = 1; });
