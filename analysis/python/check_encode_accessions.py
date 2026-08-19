import requests, json

accessions = [
    "ENCFF000BLJ",
    "ENCFF000BLK",
    "ENCFF638QZU",
    "ENCFF961FRO",
    "ENCFF001WKG",
    "ENCFF001WKD",
]

descs = {
    "ENCFF000BLJ": "K562 Repli-seq Watson strand",
    "ENCFF000BLK": "K562 Repli-seq Crick strand",
    "ENCFF638QZU": "K562 DNase-seq rep 1",
    "ENCFF961FRO": "K562 DNase-seq rep 2",
    "ENCFF001WKG": "HCT116 ChIP-seq H3K27ac",
    "ENCFF001WKD": "HCT116 ChIP-seq H3K36me3",
}

print("ENCODE REST API metadata check (/files/{ACC}/?format=json):")
print("=" * 82)

for acc in accessions:
    url = f"https://www.encodeproject.org/files/{acc}/?format=json"
    desc = descs.get(acc, "")
    try:
        r = requests.get(url, timeout=30, headers={"Accept": "application/json"})
        if r.status_code == 200:
            try:
                md = r.json()
                status = md.get("status", "?")
                ftype = md.get("file_format", "?")
                ftype_spec = md.get("file_format_type", "")
                href = md.get("href", "?")
                assembly = md.get("assembly", "?")
                replaced = md.get("replaced_by", None)
                audit = md.get("audit", {})
                flag_notes = []
                for level in ["ERROR", "WARNING", "NOT_COMPLIANT"]:
                    items = audit.get(level, [])
                    if items:
                        for it in items:
                            flag_notes.append(f"[{level}] {it.get('category','?')}: {str(it.get('detail',''))[:120]}")
                output_type = md.get("output_type", "")
                print(f"  {acc:<14} {desc:<28}")
                print(f"    HTTP 200  status=[{status}]  format={ftype}/{ftype_spec}  assembly={assembly}  output_type={output_type}")
                print(f"    href={href}")
                if replaced:
                    print(f"    ⚠️  REPLACED BY: {replaced}")
                for note in flag_notes:
                    print(f"    {note}")
                # Also test the actual href as a download
                if href.startswith("/"):
                    dl = "https://www.encodeproject.org" + href
                else:
                    dl = href
                try:
                    hr = requests.head(dl, allow_redirects=True, timeout=20)
                    print(f"    HEAD dl URL: HTTP {hr.status_code}  size={hr.headers.get('Content-Length','?')}B  type={hr.headers.get('Content-Type','')[:50]}")
                except Exception as de:
                    print(f"    HEAD dl URL ERROR: {de}")
            except Exception as e:
                print(f"  {acc:<14} {desc:<28} → HTTP 200 JSON parse error: {e}")
        elif r.status_code == 404:
            print(f"  {acc:<14} {desc:<28} → ❌ HTTP 404 (accession not found / retired)")
        else:
            print(f"  {acc:<14} {desc:<28} → HTTP {r.status_code}")
    except Exception as e:
        print(f"  {acc:<14} {desc:<28} → NET ERROR {type(e).__name__}: {e}")
    print()
