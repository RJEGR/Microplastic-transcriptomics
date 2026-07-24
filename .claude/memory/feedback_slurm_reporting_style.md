---
name: feedback-slurm-reporting-style
description: Preferencias de Ricardo para scripts SLURM y reportes técnicos — headers completos, resumen ejecutivo ≤3 frases, tablas con "Acciones recomendadas", backups antes de borrar.
metadata:
  type: feedback
---

# Estilo de scripts SLURM y reportes

## Scripts `.slurm`

**Regla:** Todo script SLURM lleva bloque de comentarios docstring en el header (propósito, entrada, salida, dependencias, dimensionado y por qué), inmediatamente después del shebang y antes de las directivas `#SBATCH`.

**Why:** Ricardo reutiliza y edita scripts meses después. Sin justificación, el "por qué" del dimensionado (CPU/RAM/time) se pierde. También sirve como paper trail cuando compara con corridas previas.

**How to apply:**
- Declarar SIEMPRE `--time`, `--mem`, `--cpus-per-task`, `--partition`, `--output`, `--error`, `--mail-user=rgomez41@uabc.edu.mx`, `--mail-type=END,FAIL`.
- Usar `set -euo pipefail` en el cuerpo.
- Rutas absolutas a `/LUSTRE/...`, nunca relativas al home.
- Arrays con `%K` explícito (concurrencia máxima) para evitar saturar el clúster compartido.
- Ver [[hpc-omica-slurm]] para el mapa partición↔etapa.

## Reportes y documentos entregables

**Regla:** Comenzar con **resumen ejecutivo de máximo 3 frases**. Tablas incluir columna **"Acciones recomendadas"**. Documentos Word: usar plantilla corporativa en `context/templates/`. Tono profesional pero accesible; evitar jerga innecesaria.

**Why:** Instrucciones globales del usuario (CLAUDE.md privado). Los stakeholders no técnicos leen sólo el ejecutivo.

**How to apply:**
- Para deliverables .md/.docx/.pptx: primer bloque = 3 frases máx.
- Tablas de hallazgos/DEGs → última columna = qué hacer con la información.
- Documentos publicables: inglés académico. Instructivos internos: español.

## Manejo de archivos

**Regla:** Nunca borrar archivo original sin crear primero copia de seguridad. Nunca compartir datos sensibles en logs.

**Why:** Instrucciones globales del usuario. En HPC compartido, un `rm` mal puesto arruina semanas de cómputo.

**How to apply:**
- Antes de sobreescribir un `.slurm`, `.R`, matriz de conteos, GTF, etc.: `cp file file.$(date +%Y%m%d).bak`.
- Rotar cualquier `NCBI_API_KEY` u otro token que aparezca en un mensaje (precedente: sesión `local_83f31324`).

Ver también [[user-profile]] y [[project-microplastics-pipeline]].
