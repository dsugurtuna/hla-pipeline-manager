FROM python:3.11-slim
LABEL maintainer="dsugurtuna"

WORKDIR /app
COPY pyproject.toml .
COPY src/ src/
RUN pip install --no-cache-dir .

ENTRYPOINT ["python", "-m", "hla_pipeline"]
