FROM dynverse/dynwrap_latest:latest

ARG GITHUB_PAT

RUN pip install scFates

RUN pip install rpy2==3.4.2 
 
RUN pip install fa2

RUN pip install git+https://github.com/dpeerlab/Palantir.git

COPY definition.yml run.py /code/

RUN ["chmod", "+x", "/code/run.py"]

ENTRYPOINT ["/code/run.py"]