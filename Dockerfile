FROM python:3.11-slim
LABEL maintainer="dsugurtuna"

RUN useradd --create-home --shell /bin/bash appuser

WORKDIR /app
COPY pyproject.toml README.md ./
COPY src/ src/
RUN pip install --no-cache-dir . \
    && chown -R appuser:appuser /app

USER appuser

# PLINK, Java and SNP2HLA are not installed here: `plan`, `verify` and
# `report` do not need them. Running a plan needs an image that has them.
ENTRYPOINT ["python", "-m", "hla_pipeline"]
CMD ["--help"]
