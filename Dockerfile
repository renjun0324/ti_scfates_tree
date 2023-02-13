FROM dynverse/dynwrap_latest:latest

ARG GITHUB_PAT

RUN pip install scFates

RUN pip install rpy2==3.4.2 
 
RUN pip install fa2

RUN pip install git+https://github.com/dpeerlab/Palantir.git

COPY definition.yml run.py example.sh /code/

RUN ["chmod", "+x", "/code/run.py"]

RUN ["chmod", "+x", "/code/example.sh"]

ENTRYPOINT ["/code/run.py"]