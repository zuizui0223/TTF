from __future__ import annotations

import json
import time
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


GBIF = "https://api.gbif.org/v1"
USER_AGENT = "ttf-relational-environment/0.2 (https://github.com/zuizui0223/TTF)"


def get_json(path: str, params: dict[str, object], retries: int = 8) -> dict:
    """Fetch one GBIF JSON endpoint with the repository's frozen retry policy."""
    url = f"{GBIF}/{path}?{urlencode(params)}"
    error: Exception | None = None
    detail = ""
    for attempt in range(retries):
        try:
            request = Request(
                url,
                headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
            )
            with urlopen(request, timeout=90) as response:
                payload = json.loads(response.read().decode("utf-8"))
            # Gentle per-request pacing. This changes no record-selection rule.
            time.sleep(0.05)
            return payload
        except HTTPError as exc:
            error = exc
            try:
                body = (
                    exc.read(512)
                    .decode("utf-8", errors="replace")
                    .replace("\n", " ")
                )
            except Exception:
                body = ""
            detail = f"HTTP {exc.code}: {body[:300]}"
            if exc.code not in {429, 500, 502, 503, 504}:
                break
            retry_after = exc.headers.get("Retry-After")
            try:
                wait = float(retry_after) if retry_after else float(2**attempt)
            except (TypeError, ValueError):
                wait = float(2**attempt)
            time.sleep(min(120.0, max(1.0, wait)))
        except Exception as exc:
            error = exc
            detail = f"{type(exc).__name__}: {exc}"
            time.sleep(min(60.0, float(2**attempt)))
    suffix = f" ({detail})" if detail else ""
    raise RuntimeError(f"GBIF request failed: {url}{suffix}") from error
