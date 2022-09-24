FROM kbase/sdkbase2:python
MAINTAINER chenry@anl.gov
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# Here we install a python coverage tool and an
# https library that is out of date in the base image.

RUN pip install --upgrade pip

# update security libraries in the base image
RUN pip install cffi --upgrade \
    && pip install pyopenssl --upgrade \
    && pip install ndg-httpsclient --upgrade \
    && pip install pyasn1 --upgrade \
    && pip install requests --upgrade \
    && pip install 'requests[security]' --upgrade \
	&& pip install --upgrade pyopenssl ndg-httpsclient && \
    pip install --upgrade pyasn1 requests 'requests[security]' && \
    pip install coverage networkx cython && \
    pip install --upgrade pip setuptools wheel cffi

RUN apt-get update
RUN apt-get install -y gcc
RUN rm -rf /miniconda/lib/python3.6/site-packages/numpy
RUN rm -rf /miniconda/lib/python3.6/site-packages/ruamel*
RUN pip install --upgrade pip
RUN pip install cobra==0.25.0
RUN pip install networkx
RUN pip install chemw==0.3.2

RUN pip install --use-deprecated=legacy-resolver git+https://github.com/ModelSEED/ModelSEEDpy.git@8a83e2ab6bb9a6ac44cff4978128a9aada1f16c6
RUN pip install git+https://github.com/Fxe/cobrakbase.git@d0967661d0056d9449621e3e1c5f01d77b6c9108

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]