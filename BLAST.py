#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VIRIONLAB — BLAST primero, GFF3 después (≥ 97% identidad)
[FIX: normalización de IDs para efetch + BANNER + RESUMEN + instalación robusta de Biopython con auto-venv]
"""

import os
import sys
import re
import platform
import shutil
import subprocess
from pathlib import Path
from datetime import datetime

# ====================== Banner ======================
VIRIONLAB_BANNER = r'''
+======================================================================+
|                                                                      |
|                           V I R I O N L A B                          |
|                                                                      |
|                         B L A S T   (pipeline)                       |
|                                                                      |
|                    Universidad Austral de Chile                      |
|                                                                      |
|                     Creado por Carlos Loncoman                       |
|                                                                      |
+======================================================================+
'''

# ====================== Parámetros por defecto ======================
DEFAULT_DB_NUCL = "nt"
DEFAULT_DB_PROT = "nr"
DEFAULT_EVALUE = "1e-20"
DEFAULT_HITLIST = 2000              # remoto (qblast) intentaremos alto; NCBI puede recortar
DEFAULT_MAX_TARGET_SEQS = "5000"    # local
DEFAULT_MIN_PIDENT = 97.0           # %
DEFAULT_MIN_QCOV = 80.0             # %
ALLOW_PROTEIN = True                # permitimos blastp si el FASTA es proteico

# ====================== Utilidades ======================

def log(msg):
    print(msg)
    LOG.append(msg)

def run(cmd, check=False):
    log(f"$ {' '.join(cmd)}")
    try:
        return subprocess.run(cmd, check=check, capture_output=True, text=True)
    except Exception as e:
        log(f"[!] Error ejecutando: {' '.join(cmd)} -> {e}")
        return None

def which(binname):
    return shutil.which(binname) is not None

# --------- AUTO-VENV para sortear PEP 668 (Homebrew / sistema gestionado) ----------
def _pip(cmd):
    """Ejecuta pip y devuelve (ok, stdout, stderr)."""
    r = subprocess.run(cmd, capture_output=True, text=True)
    return (r.returncode == 0, r.stdout, r.stderr)

def _create_and_reexec_venv(venv_dir: Path):
    """Crea venv con el mismo intérprete y re-ejecuta este script dentro del venv."""
    print(f"[i] Creando entorno virtual en: {venv_dir}")
    venv_dir.parent.mkdir(parents=True, exist_ok=True)
    # Crear venv
    r = subprocess.run([sys.executable, "-m", "venv", str(venv_dir)], capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout or "")
        print(r.stderr or "")
        print("[X] No se pudo crear el entorno virtual automáticamente.")
        return False

    # Ruta al python del venv
    vpy = venv_dir / ("Scripts/python.exe" if platform.system() == "Windows" else "bin/python")

    # Actualizar herramientas e instalar biopython
    print("[i] Instalando Biopython dentro del venv...")
    ok, out, err = _pip([str(vpy), "-m", "pip", "install", "--upgrade", "pip", "setuptools", "wheel"])
    print(out or "", end="")
    print(err or "", end="")
    ok2, out2, err2 = _pip([str(vpy), "-m", "pip", "install", "biopython>=1.85"])
    print(out2 or "", end="")
    print(err2 or "", end="")
    if not ok2:
        print("[X] No se pudo instalar Biopython dentro del venv.")
        return False

    # Re-ejecutar el script dentro del venv (evitar bucle con bandera)
    os.environ["VIRIONLAB_BOOTSTRAPPED"] = "1"
    print("[i] Relanzando el script dentro del entorno virtual...")
    os.execv(str(vpy), [str(vpy)] + sys.argv)

def ensure_biopython():
    """Garantiza que Biopython está disponible. Si el entorno está 'externally managed' (PEP 668),
    crea un venv en ~/venvs/virionlab_auto e instala allí, luego re-ejecuta el script dentro del venv.
    """
    try:
        import Bio  # noqa: F401
        return True
    except Exception:
        pass

    # Si ya relanzamos y aún no se puede importar, abortar con mensaje claro
    if os.environ.get("VIRIONLAB_BOOTSTRAPPED") == "1":
        print("[X] Biopython sigue sin estar disponible incluso dentro del venv auto-creado.")
        print("    Sugerencia: instala Xcode CLT (xcode-select --install) o usa conda/py3.11.")
        return False

    print(f"[-] Biopython no está instalado. Intentando instalar con pip en: {sys.executable}")
    print(f"[i] Plataforma: {platform.system()} {platform.machine()} | Python {platform.python_version()}")

    # 1) Herramientas básicas
    ok, out, err = _pip([sys.executable, "-m", "pip", "install", "--upgrade", "pip", "setuptools", "wheel"])
    print(out or "", end="")
    print(err or "", end="")

    # ¿Entorno gestionado externamente? (PEP 668)
    externally_managed = ("externally-managed-environment" in (err or "").lower())

    # 2) Intentar wheel binario primero
    ok2, out2, err2 = _pip([sys.executable, "-m", "pip", "install", "biopython>=1.85", "--only-binary=:all:"])
    print(out2 or "", end="")
    print(err2 or "", end="")
    if not ok2:
        # 3) Intentar instalación normal
        ok3, out3, err3 = _pip([sys.executable, "-m", "pip", "install", "biopython"])
        print(out3 or "", end="")
        print(err3 or "", end="")
        if not ok3:
            # si el error indica PEP 668 o simplemente falló, creamos venv y re-ejecutamos
            if ("externally-managed-environment" in (err3 or "").lower()) or externally_managed:
                venv_dir = Path.home() / "venvs" / "virionlab_auto"
                _create_and_reexec_venv(venv_dir)
                # si la re-ejecución falla, seguimos y devolvemos False
            else:
                print("[!] pip no pudo instalar Biopython en el entorno actual.")
                print("    Intentando crear un entorno virtual aislado...")
                venv_dir = Path.home() / "venvs" / "virionlab_auto"
                _create_and_reexec_venv(venv_dir)
                # si falla, devolver False
            return False

    # Verificar import tras instalación directa
    try:
        import Bio  # noqa: F401
        print("[+] Biopython instalado correctamente.")
        return True
    except Exception as e:
        print(f"[!] No se pudo importar Biopython tras instalar: {e}")
        # Último recurso: auto-venv
        venv_dir = Path.home() / "venvs" / "virionlab_auto"
        _create_and_reexec_venv(venv_dir)
        return False
# --------------------------------------------------------------------

def try_install_blastplus():
    osys = platform.system().lower()
    log("[-] Intentando instalar NCBI BLAST+ (mejor esfuerzo)...")
    if osys == "linux":
        if which("apt"):
            run(["sudo", "apt", "update"])
            run(["sudo", "apt", "install", "-y", "ncbi-blast+"])
            return which("blastn")
        if which("dnf"):
            run(["sudo", "dnf", "install", "-y", "ncbi-blast+"])
            return which("blastn")
        if which("yum"):
            run(["sudo", "yum", "install", "-y", "ncbi-blast+"])
            return which("blastn")
    if osys == "darwin":
        if which("brew"):
            run(["brew", "install", "blast"])
            return which("blastn")
    if osys == "windows":
        if which("choco"):
            run(["choco", "install", "ncbi-blast+", "-y"])
            return which("blastn")
    # conda/mamba
    if which("mamba"):
        run(["mamba", "install", "-y", "-c", "bioconda", "blast"])
        return which("blastn")
    if which("conda"):
        run(["conda", "install", "-y", "-c", "bioconda", "blast"])
        return which("blastn")
    log("[!] No se pudo instalar BLAST+ automáticamente.")
    return False

def desktop_outdir(base_name):
    desktop = Path.home() / "Desktop"
    safe = "".join([c for c in base_name if c.isalnum() or c in ("-", "_", " ")])
    outdir = desktop / f"resultados BLAST_{safe}_virionlab"
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir

def save_text(path, text):
    path.write_text(text, encoding="utf-8")

def append_fasta(path, header, seq):
    with open(path, "a", encoding="utf-8") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")

def write_tsv(path, header, rows):
    import csv
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for r in rows:
            w.writerow(r)

def is_nucleotide_seq(seq_str: str) -> bool:
    dna_letters = set("ACGTURYKMSWBDHVN")
    s = seq_str.upper().replace("-", "").replace(" ", "")
    if not s:
        return True
    return set(s).issubset(dna_letters)

def detect_query_type(seqs) -> str:
    flags = []
    for rec in seqs:
        flags.append("nucl" if is_nucleotide_seq(str(rec.seq)) else "prot")
    kinds = set(flags)
    if kinds == {"nucl"}:
        return "nucl"
    if kinds == {"prot"}:
        return "prot"
    return "mixed"

# === Normalización robusta de IDs para Entrez ===
ACC_RE = re.compile(r"^[A-Z]{1,4}_?\d+(\.\d+)?$", re.IGNORECASE)

def extract_accession(id_str: str) -> str:
    if not id_str:
        return id_str
    s = id_str.strip().strip("|").replace(" ", "")
    if ACC_RE.match(s):
        return s
    toks = [t for t in re.split(r"[|\s]+", s) if t]
    for t in reversed(toks):
        if ACC_RE.match(t):
            return t
    if len(toks) >= 2 and toks[1] != "gi":
        return toks[1]
    return toks[-1] if toks else s

def parse_gff3(path_gff):
    feats = {}
    with open(path_gff, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = p
            if ftype != "CDS":
                continue
            rec = {
                "seqid": seqid,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "phase": 0 if phase == "." else int(phase),
                "attrs": attrs
            }
            feats.setdefault(seqid, []).append(rec)
    return feats

def find_overlapping_cds(feats_by_seqid, sacc, sstart, send):
    if sacc not in feats_by_seqid:
        return []
    a, b = (sstart, send)
    if a > b:
        a, b = b, a
    hits = []
    for cds in feats_by_seqid[sacc]:
        c, d = cds["start"], cds["end"]
        if not (b < c or a > d):
            hits.append(cds)
    return hits

# ====================== Main ======================

if __name__ == "__main__":
    LOG = []
    print(VIRIONLAB_BANNER)
    print("Inicio:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("\n=== VIRIONLAB — BLAST primero, GFF3 después (≥97%) [FIX IDs] ===\n")

    # Biopython (robusto: import, pip o auto-venv + reexec)
    if not ensure_biopython():
        print("[X] No fue posible continuar sin Biopython.")
        sys.exit(1)

    from Bio import SeqIO
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
    from Bio import Entrez

    # BLAST+
    have_blastn = which("blastn")
    have_blastp = which("blastp")
    print(f"[i] blastn: {have_blastn} | blastp: {have_blastp}")
    if not (have_blastn or have_blastp):
        ans = input("¿Intento instalar BLAST+ automáticamente? (s/n) [s]: ").strip().lower() or "s"
        if ans in ("s", "si", "sí", "y", "yes"):
            have_blastn = try_install_blastplus() or have_blastn
            have_blastp = which("blastp") or have_blastp
            print(f"[i] blastn: {have_blastn} | blastp: {have_blastp}")

    # FASTA
    fasta_path = input("Ruta al archivo FASTA de entrada: ").strip().strip('"').strip("'")
    if not fasta_path:
        print("[X] Debes indicar un archivo FASTA.")
        sys.exit(1)
    fasta = Path(fasta_path).expanduser().resolve()
    if not fasta.exists():
        print(f"[X] No existe: {fasta}")
        sys.exit(1)

    # Carpeta salida
    outdir = desktop_outdir(fasta.stem)
    log(f"[+] Carpeta de salida: {outdir}")

    # Leer FASTA
    try:
        queries = list(SeqIO.parse(str(fasta), "fasta"))
        if not queries:
            raise ValueError("El FASTA está vacío.")
    except Exception as e:
        print(f"[X] Error leyendo FASTA: {e}")
        sys.exit(1)

    # Tipo de query
    qtype = detect_query_type(queries)
    log(f"[i] Tipo detectado: {qtype}")
    if qtype == "mixed":
        print("[X] El FASTA contiene mezcla de nucleótidos y proteínas. Usa un solo tipo.")
        sys.exit(1)
    if qtype == "prot" and not ALLOW_PROTEIN:
        print("[X] FASTA proteico detectado. Para continuar, activa ALLOW_PROTEIN=True en el script o usa FASTA nucleótido.")
        sys.exit(1)

    # Parámetros de filtro y búsqueda
    try:
        min_pident = float(input(f"Umbral identidad mínima [%] [{DEFAULT_MIN_PIDENT}]: ").strip() or DEFAULT_MIN_PIDENT)
    except:
        min_pident = DEFAULT_MIN_PIDENT
    try:
        min_qcov = float(input(f"Umbral cobertura mínima de la query [%] [{DEFAULT_MIN_QCOV}]: ").strip() or DEFAULT_MIN_QCOV)
    except:
        min_qcov = DEFAULT_MIN_QCOV

    mode = "remoto"
    if (qtype == "nucl" and have_blastn) or (qtype == "prot" and have_blastp):
        sel = input("Modo de búsqueda: [1] Remoto NCBI  [2] Local\nElige 1/2 [1]: ").strip() or "1"
        mode = "local" if sel == "2" else "remoto"
    else:
        print("[i] No hay ejecutable local apropiado; se usará BLAST remoto.")
        mode = "remoto"

    # DB y programa
    if qtype == "nucl":
        program = "blastn"
        db_remote = DEFAULT_DB_NUCL
        blast_bin = "blastn"
    else:
        program = "blastp"
        db_remote = DEFAULT_DB_PROT
        blast_bin = "blastp"

    # Entrez (para efetch)
    email = input("Email para Entrez (obligatorio para descargar FASTA): ").strip()
    if not email:
        print("[X] Debes proporcionar un email para Entrez.")
        sys.exit(1)
    api_key = input("API key de NCBI (opcional, Enter para omitir): ").strip() or None
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    # Rutas de salida
    tsv_all = outdir / "blast_all.tsv"
    tsv_filt = outdir / "blast_filtered.tsv"
    fasta_hits = outdir / "hits_>=97pc.fasta"
    annot_csv = outdir / "annotations_from_gff3.tsv"
    blast_xml = outdir / "blast.xml"   # si remoto, guardamos XML
    for p in (tsv_all, tsv_filt, fasta_hits, annot_csv, blast_xml):
        if p.exists():
            p.unlink()

    # Campos para outfmt 6
    header = ["qseqid","sacc","stitle","pident","length","evalue","bitscore","qstart","qend","sstart","send","qcovhsp"]
    outfmt6_fields = "6 " + " ".join(header)

    # === Ejecutar BLAST ===
    rows_all = []
    if mode == "local":
        # Guardar queries a fichero temporal
        qf = outdir / "queries_tmp.fasta"
        if qf.exists():
            qf.unlink()
        for rec in queries:
            append_fasta(qf, rec.id, str(rec.seq))

        # Ejecutar local (tabular)
        db_pref = input("Prefijo de base local (resultado de makeblastdb): ").strip()
        cmd = [
            blast_bin,
            "-query", str(qf),
            "-db", db_pref,
            "-outfmt", outfmt6_fields,
            "-evalue", DEFAULT_EVALUE,
            "-max_target_seqs", DEFAULT_MAX_TARGET_SEQS,
            "-num_threads", str(os.cpu_count() or 1),
        ]
        if program == "blastn":
            cmd += ["-task", "megablast"]
            try:
                pid_arg = input(f"¿Aplicar -perc_identity {int(min_pident)} en BLAST local? (s/n) [s]: ").strip().lower() or "s"
            except:
                pid_arg = "s"
            if pid_arg.startswith("s"):
                cmd += ["-perc_identity", str(int(min_pident))]

        r = run(cmd)
        if r and r.stdout:
            for line in r.stdout.strip().splitlines():
                parts = line.split("\t")
                if len(parts) >= len(header):
                    parts[1] = extract_accession(parts[1])  # NORMALIZAR sacc
                    rows_all.append(parts)
        else:
            # fallback: a archivo
            cmd_out = cmd + ["-out", str(tsv_all)]
            r2 = run(cmd_out)
            if not r2 or r2.returncode != 0:
                print("[X] BLAST local falló.")
                save_text(outdir / "log.txt", "\n".join(LOG))
                sys.exit(1)
            for ln in tsv_all.read_text().splitlines():
                if not ln.strip():
                    continue
                parts = ln.split("\t")
                if len(parts) >= len(header):
                    parts[1] = extract_accession(parts[1])  # NORMALIZAR sacc
                    rows_all.append(parts)
    else:
        # Remoto qblast -> XML, parsear a filas tipo outfmt 6
        log(f"[i] BLAST remoto: program={program}, db={db_remote}, hitlist={DEFAULT_HITLIST}")
        xml_strings = []
        for rec in queries:
            log(f"[.] qblast -> {rec.id} (len={len(rec.seq)})")
            try:
                handle = NCBIWWW.qblast(
                    program=program,
                    database=db_remote,
                    sequence=str(rec.seq),
                    expect=float(DEFAULT_EVALUE),
                    hitlist_size=DEFAULT_HITLIST,
                    megablast=(program == "blastn")
                )
                xml = handle.read()
                handle.close()
                xml_strings.append(xml)
            except Exception as e:
                log(f"[!] qblast error en {rec.id}: {e}")

        if not xml_strings:
            print("[X] No se obtuvieron resultados desde NCBI.")
            save_text(outdir / "log.txt", "\n".join(LOG))
            sys.exit(1)

        save_text(blast_xml, "\n".join(xml_strings))
        log("[+] XML remoto guardado. Parseando...")

        with open(blast_xml, "r", encoding="utf-8") as fh:
            records = NCBIXML.parse(fh)
            for record in records:
                qid = getattr(record, "query_id", None) or getattr(record, "query", "query_unk").split()[0]
                qlen = getattr(record, "query_length", None)
                for aln in record.alignments:
                    sacc_raw = getattr(aln, "accession", None) or aln.hit_id
                    sacc = extract_accession(sacc_raw)  # NORMALIZAR
                    stitle = getattr(aln, "hit_def", "")
                    for hsp in aln.hsps:
                        ql = qlen
                        if not ql:
                            qrec = next((x for x in queries if x.id == qid), None)
                            ql = len(qrec.seq) if qrec else None
                        qcov = round((hsp.align_length / ql) * 100.0, 3) if (ql and hsp.align_length) else ""
                        pident = round((hsp.identities / hsp.align_length) * 100.0, 3) if hsp.align_length else ""
                        row = [
                            str(qid),
                            str(sacc),
                            str(stitle),
                            str(pident),
                            str(hsp.align_length),
                            str(hsp.expect),
                            str(hsp.bits),
                            str(hsp.query_start),
                            str(hsp.query_end),
                            str(hsp.sbjct_start),
                            str(hsp.sbjct_end),
                            str(qcov)
                        ]
                        rows_all.append(row)

    # Guardar todo (sin filtrar)
    write_tsv(tsv_all, header, rows_all)
    log(f"[+] BLAST bruto guardado: {tsv_all}")

    # === Filtrado por identidad y cobertura ===
    rows_f = []
    for r in rows_all:
        try:
            pid = float(r[3])
            qcov = float(r[11])
        except:
            continue
        if pid >= min_pident and qcov >= min_qcov:
            rows_f.append(r)
    write_tsv(tsv_filt, header, rows_f)
    log(f"[+] Filtrados (pident>={min_pident} & qcov>={min_qcov}): {len(rows_f)} filas -> {tsv_filt}")

    # Si no hay filtrados, cerramos con resumen básico
    if not rows_f:
        print("\n— Resumen —")
        print(f"  • Queries procesadas: {len(queries)}")
        print(f"  • Hits totales (sin filtrar): {len(rows_all)}")
        print(f"  • Hits filtrados (≥{min_pident}% identidad y ≥{min_qcov}% qcovhsp): 0")
        print("\nGracias por usar este script de VIRIONLAB 🙌")
        print("Sugerencia: relaja los umbrales (p.ej., identidad 95% o qcovhsp 60–70%) y reintenta.")
        save_text(outdir / "log.txt", "\n".join(LOG))
        sys.exit(0)

    # === Descargar FASTA de sujetos que pasaron ===
    # sacc únicos y NORMALIZADOS
    saccs = []
    seen = set()
    for r in rows_f:
        acc = extract_accession(r[1])
        if acc and acc not in seen:
            seen.add(acc)
            saccs.append(acc)

    log(f"[i] Descargando FASTA de {len(saccs)} accesiones (en lotes).")
    db_entrez = "nucleotide" if qtype == "nucl" else "protein"

    def batched(lst, n=200):
        for i in range(0, len(lst), n):
            yield lst[i:i+n]

    if fasta_hits.exists():
        fasta_hits.unlink()

    for batch in batched(saccs, 200):
        ids = ",".join(batch)  # SIN pipes, solo accession.version
        try:
            with Entrez.efetch(db=db_entrez, id=ids, rettype="fasta", retmode="text") as h:
                fasta_text = h.read()
            if fasta_text.strip():
                with open(fasta_hits, "a", encoding="utf-8") as outfa:
                    outfa.write(fasta_text if fasta_text.endswith("\n") else fasta_text + "\n")
        except Exception as e:
            log(f"[!] efetch error en batch (n={len(batch)}): {e}")

    log(f"[+] FASTA de hits guardado: {fasta_hits}")

    # === (Opcional) Anotar con GFF3 (overlap HSP ↔ CDS) ===
    use_gff = input("¿Deseas anotar HSPs con un GFF3 (sacc==seqid)? (s/n) [n]: ").strip().lower() or "n"
    if use_gff.startswith("s"):
        gff_path = input("Ruta al archivo GFF3: ").strip().strip('"').strip("'")
        gff = Path(gff_path).expanduser().resolve()
        if not gff.exists():
            log("[!] GFF3 no encontrado. Omitiendo anotación.")
        else:
            feats_by_seqid = parse_gff3(str(gff))
            ann_rows = []
            ann_header = header + ["cds_seqid","cds_start","cds_end","cds_strand","cds_phase","cds_attrs"]
            for r in rows_f:
                sacc = extract_accession(r[1])  # por si el GFF usa accession limpio
                try:
                    sstart = int(r[9]); send = int(r[10])
                except:
                    continue
                overlaps = find_overlapping_cds(feats_by_seqid, sacc, sstart, send)
                if not overlaps:
                    ann_rows.append(r + ["", "", "", "", "", ""])
                else:
                    for cds in overlaps:
                        ann_rows.append(r + [
                            cds["seqid"], cds["start"], cds["end"], cds["strand"], cds["phase"], cds["attrs"]
                        ])
            write_tsv(annot_csv, ann_header, ann_rows)
            log(f"[+] Anotación GFF3 guardada: {annot_csv}")

     # === Resumen final ===
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    unique_hits = len(set(extract_accession(r[1]) for r in rows_f))

    print("\n— Resumen —")
    print(f"  • Queries procesadas: {len(queries)}")
    print(f"  • Hits totales (sin filtrar): {len(rows_all)}")
    print(f"  • Hits filtrados: {len(rows_f)}")
    print(f"    ↳ Criterio de filtración final:")
    print(f"       - Identidad (pident) ≥ {min_pident}%")
    print(f"       - Cobertura de la query (qcovhsp) ≥ {min_qcov}%")
    print("         * qcovhsp se calcula como (longitud del HSP / longitud de la query) × 100.")
    if mode == "remoto":
        print("         * Nota: en modo remoto (qblast/XML) qcovhsp es una aproximación derivada de hsp.align_length y query_length.")

    print(f"\n  • Accesiones únicas descargadas (FASTA): {unique_hits}")
    print("    ↳ Antes de la descarga FASTA se deduplican accesiones (sólo se baja una vez cada 'sacc').")

    print("\n  • Archivos generados:")
    print(f"      - {tsv_all.name}       (todos los HSPs sin filtrar)")
    print(f"      - {tsv_filt.name}      (sólo HSPs que cumplen pident y qcovhsp)")
    print(f"      - {fasta_hits.name}    (secuencias FASTA de accesiones únicas filtradas)")
    if (outdir / 'annotations_from_gff3.tsv').exists():
        print(f"      - {annot_csv.name}     (anotación de solapamiento HSP ↔ CDS desde GFF3, si se proporcionó)")
    if blast_xml.exists():
        print(f"      - {blast_xml.name}     (salida XML de BLAST en modo remoto)")

    print("\n— ¿Por qué 'hits totales' puede ser mayor que 'hits filtrados'? —")
    print("  • Incluyen HSPs con identidad o cobertura por debajo del umbral.")
    print("  • Puede haber múltiples HSPs por el mismo sujeto; sólo se retienen los que cumplen los umbrales.")
    print("  • La descarga FASTA usa accesiones únicas; esto reduce conteos respecto al total de HSPs.")

    print("\n— ¿Qué significan los umbrales? —")
    print(f"  • Identidad ≥ {min_pident}%: de cada 100 posiciones alineadas, al menos {min_pident} son idénticas.")
    print(f"  • Cobertura de la query (qcovhsp) ≥ {min_qcov}%: el HSP cubre al menos el {min_qcov}% de la longitud de tu secuencia de consulta.")
    print("    - Ej.: si tu query mide 1000 nt y qcovhsp=80%, el alineamiento abarca ~800 nt de tu query.")

    print("\nGracias por usar este script de VIRIONLAB 🙌")
    print(f"Fin: {timestamp}")

    # Guardar log
    save_text(outdir / "log.txt", "\n".join(LOG))

    print("\n✅ Revisa tu Escritorio en la carpeta:")
    print(f"   {outdir}")
