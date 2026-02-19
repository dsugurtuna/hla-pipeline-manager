.PHONY: install dev test lint clean docker

install:
	pip install -e .

dev:
	pip install -e ".[dev]"

test:
	pytest tests/ -v

lint:
	ruff check src/ tests/

clean:
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	rm -rf .pytest_cache dist *.egg-info

docker:
	docker build -t hla-pipeline-manager .
