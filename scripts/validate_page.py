#!/usr/bin/env python3
"""End-to-end browser checks for the static lookup page."""

from __future__ import annotations

from pathlib import Path

from playwright.sync_api import expect, sync_playwright


URL = "http://127.0.0.1:4173/"
SCREENSHOT = Path("/tmp/lookup-page-validation.png")


def main() -> None:
    page_errors: list[str] = []
    console_errors: list[str] = []

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        context = browser.new_context(viewport={"width": 1440, "height": 1000})
        context.grant_permissions(
            ["clipboard-read", "clipboard-write"],
            origin=URL.rstrip("/"),
        )
        page = context.new_page()
        page.on("pageerror", lambda exc: page_errors.append(str(exc)))
        page.on(
            "console",
            lambda msg: console_errors.append(msg.text) if msg.type == "error" else None,
        )

        page.goto(URL, wait_until="domcontentloaded", timeout=30000)
        expect(page.locator("#datasetMeta")).to_contain_text(
            "2,464 EC buckets",
            timeout=10000,
        )
        expect(page.locator("#rdkitStatus")).to_contain_text(
            "RDKit.js",
            timeout=60000,
        )

        page.locator(".browse-item").filter(has_text="Oxidoreductases").first.click()
        expect(page.locator("#breadcrumbs")).to_contain_text("1")
        page.locator(".browse-item").filter(has_text="Acting on the CH-OH group").first.click()
        expect(page.locator("#breadcrumbs")).to_contain_text("1.1")
        page.locator(".browse-item").filter(
            has_text="With NAD(+) or NADP(+) as acceptor",
        ).first.click()
        expect(page.locator("#breadcrumbs")).to_contain_text("1.1.1")
        page.locator(".browse-item").filter(has_text="1.1.1.1").first.click()
        expect(page.locator("#panelHead")).to_contain_text(
            "alcohol dehydrogenase",
            timeout=10000,
        )
        expect(page.locator(".template-card").first).to_be_visible(timeout=10000)
        expect(page.locator(".reaction-art svg").first).to_be_visible(timeout=60000)
        expect(page.locator(".smarts-block pre").first).to_contain_text(">>")

        first_copy = page.locator(".copy-template").first
        first_smarts = first_copy.get_attribute("data-copy")
        first_copy.click()
        expect(first_copy).to_have_text("Copied", timeout=5000)
        copied = page.evaluate("navigator.clipboard.readText()")
        assert copied == first_smarts, "clipboard SMARTS did not match first template"

        page.locator("#searchInput").fill("[O&H1&+0&D1:1]>>")
        expect(page.locator("#panelHead")).to_contain_text(
            "template matches",
            timeout=10000,
        )
        assert page.locator(".template-card").count() > 0, "SMARTS search returned no templates"
        expect(page.locator(".reaction-art svg").first).to_be_visible(timeout=60000)

        page.screenshot(path=str(SCREENSHOT), full_page=True)
        browser.close()

    assert not page_errors, page_errors
    bad_console = [msg for msg in console_errors if "Failed to load resource" not in msg]
    assert not bad_console, bad_console[:5]
    print("browser ok")
    print(f"screenshot {SCREENSHOT}")


if __name__ == "__main__":
    main()
