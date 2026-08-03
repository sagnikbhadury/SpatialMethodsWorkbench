const { chromium } = require('C:/tmp/smw-video-tools/node_modules/playwright-core');
const path = require('path');

const outputDir = path.resolve('media/youtube/captures');

async function capture(page, name) {
  await page.waitForTimeout(1200);
  await page.screenshot({ path: path.join(outputDir, name), fullPage: false });
}

async function openTab(page, name) {
  await page.getByRole('tab', { name, exact: true }).click();
  await page.waitForTimeout(800);
}

(async () => {
  const browser = await chromium.launch({
    executablePath: 'C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe',
    headless: true
  });
  const page = await browser.newPage({ viewport: { width: 1920, height: 1080 }, deviceScaleFactor: 1 });
  await page.goto('https://sagnikbhadury.shinyapps.io/spatial-methods-workbench/', { waitUntil: 'domcontentloaded', timeout: 120000 });
  await page.waitForSelector('#run', { timeout: 120000 });
  await page.waitForTimeout(5000);

  await capture(page, '01_home_and_mapping.png');

  await openTab(page, 'Readiness');
  await page.locator('#readiness').scrollIntoViewIfNeeded();
  await capture(page, '02_readiness.png');

  await openTab(page, 'Usage & restrictions');
  await page.getByRole('heading', { name: 'Restrictions and responsibilities' }).scrollIntoViewIfNeeded();
  await capture(page, '03_restrictions.png');

  await openTab(page, 'Pipeline reference');
  await page.getByRole('heading', { name: 'Common programmatic interface' }).scrollIntoViewIfNeeded();
  await capture(page, '04_pipeline_reference.png');

  await page.locator('#run').scrollIntoViewIfNeeded();
  await capture(page, '05_controls_and_citation.png');

  await openTab(page, 'Usage & restrictions');
  await page.getByRole('heading', { name: 'Install the complete suite locally' }).scrollIntoViewIfNeeded();
  await capture(page, '09_installation.png');

  await openTab(page, 'Discuss results');
  await capture(page, '10_collaboration_boundary.png');

  await browser.close();
})();
