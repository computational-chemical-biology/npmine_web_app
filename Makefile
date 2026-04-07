ENV_FILE ?= .env
COMPOSE := docker compose --env-file $(ENV_FILE) -f docker-compose.yml
FLASK := $(COMPOSE) exec web conda run -n npmine_web_app flask
SHELL := /bin/bash

.PHONY: help up down restart logs ps shell migrate ensure-web import-compounds import-compounds-enriched import-taxa

help:
	@echo "Available targets:"
	@echo "  make up                         Start the app stack"
	@echo "  make down                       Stop the app stack"
	@echo "  make restart                    Restart the app stack"
	@echo "  make logs                       Follow compose logs"
	@echo "  make ps                         Show running services"
	@echo "  make shell                      Open a shell in the web container"
	@echo "  make migrate                    Run database migrations"
	@echo "  make import-compounds           Import compounds without external enrichment"
	@echo "  make import-compounds-enriched  Import compounds with PubChem and NPClassifier"
	@echo "  make import-taxa                Import taxa"

up:
	$(COMPOSE) up --build -d

down:
	$(COMPOSE) down

restart: down up

logs:
	$(COMPOSE) logs -f

ps:
	$(COMPOSE) ps

shell:
	$(COMPOSE) exec web bash

migrate:
	$(FLASK) db upgrade

ensure-web:
	$(COMPOSE) up --build -d

import-compounds: ensure-web
	$(FLASK) import-compounds

import-compounds-enriched: ensure-web
	$(FLASK) import-compounds --with-pubchem --with-npclassifier

import-taxa: ensure-web
	$(FLASK) import-taxa
