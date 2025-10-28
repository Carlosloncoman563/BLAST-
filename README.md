# BLAST-
Pipeline automatizado de VIRIONLAB para ejecutar BLAST (local o remoto), filtrar resultados por ≥97 % de identidad y ≥80 % de cobertura, descargar secuencias FASTA desde NCBI y, opcionalmente, anotar solapamientos con un archivo GFF3. Incluye instalación automática de Biopython y generación de carpeta de salida en el Escritorio.
Este script implementa un pipeline robusto denominado “VIRIONLAB — BLAST primero, GFF3 después (≥97% identidad)” para buscar homologías de secuencias y, opcionalmente, anotar los alineamientos significativos contra características codificantes (CDS) definidas en un archivo GFF3. Está pensado para uso práctico en laboratorio (diagnóstico/investigación) y prioriza la reproducibilidad, la tolerancia a errores y la portabilidad.
Al iniciar, muestra un banner informativo, registra un log de todas las acciones y valida la disponibilidad de Biopython. Si no está instalado o si el entorno es “externally managed” (PEP 668), el script crea automáticamente un entorno virtual en ~/venvs/virionlab_auto, instala Biopython y re-ejecuta el propio script dentro de ese entorno, evitando interrupciones. Adicionalmente, detecta la presencia de BLAST+ (blastn/blastp) y, si falta, intenta instalarlo con los gestores disponibles (apt/dnf/yum, Homebrew, choco o conda/mamba), sin bloquear la ejecución si no se logra.
El usuario proporciona un FASTA de entrada. El script detecta el tipo de secuencia (nucleótido/proteína), permite usar BLAST remoto (NCBI qblast) o local (si hay base creada con makeblastdb), y aplica parámetros por defecto razonables (e.g., e-value 1e-20, megablast para nucleótidos). Los resultados brutos se capturan en formato tabular tipo outfmt 6 (o se parsean desde XML en modo remoto) con campos clave (qseqid, sacc, stitle, pident, qcovhsp, coordenadas, etc.). El script normaliza accesiones (sacc) para uso consistente con Entrez (efetch), filtra por identidad mínima (≥ 97 %) y cobertura de la query (≥ 80 %, ambos configurables), y genera salidas claras:
blast_all.tsv (todos los HSPs),
blast_filtered.tsv (solo hits que cumplen umbrales),
hits_>=97pc.fasta (descargado por Entrez a partir de accesiones únicas),
blast.xml (si se usó qblast),
annotations_from_gff3.tsv (opcional).
Si el usuario lo desea, el script carga un GFF3 y calcula solapamientos HSP↔CDS por seqid y coordenadas, añadiendo columnas con strand, phase y attrs. Todos los artefactos se guardan en una carpeta auto-creada en el Escritorio: resultados BLAST_<basename>_virionlab. Al finalizar, imprime un resumen interpretativo (conteos totales vs. filtrados, criterios aplicados, notas sobre qcovhsp) y deja un log con los comandos ejecutados, facilitando auditoría, trazabilidad y repetición del análisis.
VIRIONLAB — BLAST primero, GFF3 después (≥97% identidad)
Pipeline de línea de comando para buscar homologías con BLAST (remoto NCBI o local con BLAST+), filtrar por identidad y cobertura, descargar FASTA desde NCBI (Entrez/efetch) y anotar solapamientos HSP↔CDS usando un GFF3. Incluye auto-instalación de Biopython con auto-venv cuando el entorno está “externally managed” (PEP 668).
Salida organizada automáticamente en el Escritorio: resultados BLAST_<basename>_virionlab/.
👀 Qué hace
Valida/instala Biopython. Si el entorno es gestionado (PEP 668), crea un venv en ~/venvs/virionlab_auto e inicia el script allí.
(Opcional) Instala BLAST+ con apt/dnf/yum/Homebrew/choco/conda/mamba (mejor esfuerzo).
Lee un FASTA, detecta si es nucleótido o proteína.
Ejecuta BLAST remoto (qblast) o local (si hay blastn/blastp y DB local).
Normaliza accesiones (sacc) para Entrez.
Filtra por pident ≥ 97 % y qcovhsp ≥ 80 % (configurable en tiempo de ejecución).
Descarga FASTA de accesiones únicas mediante Entrez/efetch.
(Opcional) Anota con GFF3 los solapamientos HSP↔CDS.
Genera TSV/FASTA/XML y un log reproducible.
🧩 Requisitos
Python 3.9+ (recomendado 3.11/3.12).
Internet (para qblast/Entrez).
Para BLAST local: NCBI BLAST+ y una base creada con makeblastdb.
Permisos de instalación si el script intenta instalar BLAST+.
El script se encarga de Biopython. No necesitas instalarlo manualmente.
📦 Archivos de salida
En ~/Desktop/resultados BLAST_<basename>_virionlab/:
blast_all.tsv → HSPs sin filtrar (formato tipo outfmt 6).
blast_filtered.tsv → sólo HSPs que pasan los umbrales.
hits_>=97pc.fasta → secuencias de sujetos filtrados (accesiones únicas).
blast.xml → salida cruda de qblast (si usaste modo remoto).
annotations_from_gff3.tsv → solapamientos HSP↔CDS (si proporcionas GFF3).
log.txt → comandos/acciones ejecutados.
▶️ Ejecución rápida (interactiva)
python virionlab_blast_gff3.py
El script te irá preguntando:
Ruta al FASTA
Umbrales (identidad % y qcov %)
Modo: Remoto (NCBI) o Local (si hay blastn/blastp)
Email Entrez (obligatorio) y API key (opcional, acelera cuotas)
(Local) Prefijo de DB (-db) creado con makeblastdb
(Opcional) Ruta a GFF3 para anotar
Al final verás un resumen y la ruta de la carpeta de resultados.
🧪 Ejemplo completo (BLAST remoto)
Supón que tienes prv_s1.fasta (nucleótidos):
Ruta al archivo FASTA de entrada: /Users/carlos/Desktop/prv_s1.fasta
Umbral identidad mínima [%] [97.0]: 97
Umbral cobertura mínima de la query [%] [80.0]: 80
Modo de búsqueda: [1] Remoto NCBI  [2] Local
Elige 1/2 [1]: 1
Email para Entrez (obligatorio para descargar FASTA): tu@uach.cl
API key de NCBI (opcional, Enter para omitir): (Enter)
¿Deseas anotar HSPs con un GFF3 (sacc==seqid)? (s/n) [n]: s
Ruta al archivo GFF3: /Users/carlos/ref/PRV.gff3
Resultado:
resultados BLAST_prv_s1_virionlab/ con blast_all.tsv, blast_filtered.tsv, hits_>=97pc.fasta, blast.xml, annotations_from_gff3.tsv, log.txt.
🧰 Ejemplo (BLAST local)
Crea la DB (una sola vez):
makeblastdb -in ref_prv_nt.fasta -dbtype nucl -out db/prv_nt
Ejecuta el script y elige Local:
Ruta al archivo FASTA de entrada: /Users/carlos/Desktop/prv_queries.fasta
...
Elige 1/2 [1]: 2
Prefijo de base local (resultado de makeblastdb): db/prv_nt
¿Aplicar -perc_identity 97 en BLAST local? (s/n) [s]: s
El script usará -task megablast si es nucleótido y añadirá -perc_identity si aceptas.
🎛️ Parámetros y comportamiento
Umbrales por defecto: pident=97.0, qcovhsp=80.0. Puedes ajustarlos cuando te los pregunte.
E-value por defecto: 1e-20.
Hitlist remoto: 2000 (NCBI puede recortar).
Máximo targets local: 5000 (-max_target_seqs).
Normalización de accesiones: limpia variantes de IDs compuestos para que Entrez acepte acc.version.
qcovhsp (remoto): se aproxima como (hsp.align_length / query_length) × 100.
GFF3: se asume que seqid corresponde al sacc (accesión); se listan todas las CDS que solapan con el rango del HSP.
🛠️ Consejos y buenas prácticas
Email Entrez real y API key (si la tienes) mejoran estabilidad y cupos de NCBI.
Si el FASTA mezcla nucleótido y proteína, el script aborta: usa archivos separados.
En local, verifica que la DB y el tipo coincidan (nucl↔blastn, prot↔blastp).
Si fallan instalaciones del sistema (apt/brew/choco), instala BLAST+ manualmente o via conda install -c bioconda blast.
Guarda el GFF3 con seqid igual a accesión (p. ej., NC_XXXX.1) si quieres anotación exacta.
🐞 Solución de problemas
“Biopython no disponible”: el script creará ~/venvs/virionlab_auto y se re-ejecutará. Si aún falla, instala Xcode CLT (macOS) o usa conda.
“No hay blastn/blastp”: elige Remoto o instala BLAST+ (ver arriba).
“Cero filtrados”: relaja umbrales (ej. identidad 95 % o qcovhsp 60–70 %).
“GFF3 sin coincidencias”: revisa que el seqid del GFF3 coincida con sacc (misma accesión y versión).
📄 Licencia y crédito
Desarrollado por VirionLab (UACh) — preparado para rutinas de diagnóstico e investigación en virología.
Incluye registro de comandos para trazabilidad y reproducibilidad.
