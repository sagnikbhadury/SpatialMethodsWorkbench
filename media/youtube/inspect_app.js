const { chromium } = require('C:/tmp/smw-video-tools/node_modules/playwright-core');

(async () => {
  const browser = await chromium.launch({
    executablePath: 'C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe',
    headless: true
  });
  const page = await browser.newPage({ viewport: { width: 1920, height: 1080 }, deviceScaleFactor: 1 });
  await page.goto('https://sagnikbhadury.shinyapps.io/spatial-methods-workbench/', { waitUntil: 'domcontentloaded', timeout: 120000 });
  await page.waitForSelector('#run', { timeout: 120000 });
  await page.waitForTimeout(5000);
  const info = await page.evaluate(() => ({
    title: document.title,
    tabs: [...document.querySelectorAll('[role="tab"], .nav-link')].map(x => ({ text: x.innerText.trim(), href: x.getAttribute('href') })),
    inputs: [...document.querySelectorAll('input,select,button')].map(x => ({ tag: x.tagName, id: x.id, type: x.type, text: (x.innerText || '').trim(), value: x.value })).filter(x => x.id || x.text),
    headings: [...document.querySelectorAll('h1,h2,h3')].map(x => x.innerText.trim())
  }));
  process.stdout.write(JSON.stringify(info, null, 2));
  await page.screenshot({ path: 'media/youtube/initial.png', fullPage: false });
  await browser.close();
})();
