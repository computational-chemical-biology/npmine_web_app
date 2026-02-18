FROM continuumio/miniconda3:latest AS builder

WORKDIR /app
COPY environment.yml environment.yml
RUN conda env create -f environment.yml --name npmine_web_app

FROM continuumio/miniconda3:latest AS production

WORKDIR /app

COPY --from=builder /opt/conda/envs/npmine_web_app /opt/conda/envs/npmine_web_app
COPY . .

ENV PATH=/opt/conda/envs/npmine_web_app/bin:$PATH

EXPOSE 5000

CMD ["gunicorn", "--bind", "0.0.0.0:5000", "websiteNPMINE:create_app()"]
