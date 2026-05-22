#!/usr/bin/env python3

import os
import pandas as pd
import subprocess
import zipfile
import shutil

INPUT_FILE = "test.csv"   # cambia esto por tu archivo
OUTPUT_DIR = "OUTPUT"     # carpeta raíz de salida
FILETYPE = "csv"          # "csv" o "tsv"

FAILED_REPORT = os.path.join(OUTPUT_DIR, "failed_downloads.tsv")


def run_cmd(cmd):
    print(" ".join(cmd))
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip())
    return result.stdout


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    if FILETYPE == "csv":
        df = pd.read_csv(INPUT_FILE)
    elif FILETYPE == "tsv":
        df = pd.read_csv(INPUT_FILE, sep="\t")
    else:
        raise ValueError("FILETYPE must be 'csv' or 'tsv'")

    if "HOG" not in df.columns or "Protein" not in df.columns:
        raise ValueError("El archivo debe tener columnas llamadas exactamente: HOG y Protein")

    failed = []

    for i, row in df.iterrows():
        hog = str(row["HOG"]).strip()
        prot = str(row["Protein"]).strip()

        if hog == "nan" or prot == "nan":
            failed.append((hog, prot, "Missing HOG or Protein value"))
            continue

        outdir = os.path.join(OUTPUT_DIR, hog, prot, "seqs")
        os.makedirs(outdir, exist_ok=True)

        zip_path = os.path.join(outdir, f"{prot}.zip")
        tmp_extract = os.path.join(outdir, "tmp_extract")

        protein_out = os.path.join(outdir, "protein.faa")
        rna_out = os.path.join(outdir, "rna.fna")

        # Si ya existen ambos, skip
        if os.path.exists(protein_out) and os.path.exists(rna_out):
            print(f"[SKIP] {prot} ya descargada")
            continue

        print(f"[DOWNLOADING] HOG={hog} Protein={prot}")

        cmd = [
            "datasets", "download", "gene", "accession", prot,
            "--filename", zip_path
        ]

        try:
            run_cmd(cmd)
        except Exception as e:
            failed.append((hog, prot, f"Download failed: {e}"))
            print(f"[ERROR] Descarga fallida {prot}")
            continue

        # Extraer zip
        try:
            if os.path.exists(tmp_extract):
                shutil.rmtree(tmp_extract)
            os.makedirs(tmp_extract)

            with zipfile.ZipFile(zip_path, 'r') as z:
                z.extractall(tmp_extract)

        except Exception as e:
            failed.append((hog, prot, f"Extraction failed: {e}"))
            print(f"[ERROR] Extracción fallida {prot}")
            continue

        # Buscar archivos protein y rna dentro del zip
        protein_file = None
        rna_file = None

        for root, dirs, files in os.walk(tmp_extract):
            for f in files:
                fname = f.lower()

                # típicamente: protein.faa / rna.fna
                if fname.startswith("protein") and fname.endswith((".faa", ".fasta", ".fa")):
                    protein_file = os.path.join(root, f)

                if fname.startswith("rna") and fname.endswith((".fna", ".fasta", ".fa")):
                    rna_file = os.path.join(root, f)

        if protein_file is None:
            failed.append((hog, prot, "Protein file not found inside ZIP"))
            print(f"[ERROR] No se encontró archivo protein para {prot}")
            shutil.rmtree(tmp_extract)
            continue

        if rna_file is None:
            failed.append((hog, prot, "RNA file not found inside ZIP"))
            print(f"[ERROR] No se encontró archivo rna para {prot}")
            shutil.rmtree(tmp_extract)
            continue

        # Copiar ambos ficheros al directorio final
        try:
            shutil.copy(protein_file, protein_out)
            shutil.copy(rna_file, rna_out)
        except Exception as e:
            failed.append((hog, prot, f"Copy failed: {e}"))
            print(f"[ERROR] No se pudieron copiar archivos protein/rna para {prot}")
            shutil.rmtree(tmp_extract)
            continue

        shutil.rmtree(tmp_extract)

        print(f"[OK] Guardado protein en {protein_out}")
        print(f"[OK] Guardado rna en {rna_out}")

    # Guardar reporte final
    if failed:
        failed_df = pd.DataFrame(failed, columns=["HOG", "Protein", "Reason"])
        failed_df.to_csv(FAILED_REPORT, sep="\t", index=False)
        print(f"\n[REPORT] Proteínas fallidas guardadas en: {FAILED_REPORT}")
    else:
        print("\n[REPORT] Todas las proteínas se descargaron correctamente. No hay fallos.")


if __name__ == "__main__":
    main()