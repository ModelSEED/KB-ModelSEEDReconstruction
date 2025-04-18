FROM kbase/sdkpython:3.8.0
MAINTAINER chenry@anl.gov
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# Here we install a python coverage tool and an
# https library that is out of date in the base image.

RUN pip install --upgrade pip
RUN pip install chemw==0.3.2
RUN pip install pandas

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
RUN pip install "numpy<1.24"
RUN pip install cobra
RUN pip install networkx
RUN pip install deepdiff
RUN pip install h5py

RUN echo '38' >/dev/null && pip install --use-deprecated=legacy-resolver git+https://github.com/cshenry/ModelSEEDpy.git

RUN echo '31' >/dev/null && mkdir deps && cd deps && \
	git clone --branch main https://github.com/cshenry/KBBaseModules.git
RUN echo '1' >/dev/null && cd deps && \
	git clone https://github.com/ModelSEED/ModelSEEDDatabase.git && \
    cd ModelSEEDDatabase && git checkout 3346b71a34bc9d8c5a365b71d5a2959ffbe6c26e
RUN echo '2' >/dev/null && cd deps && \
    git clone https://github.com/cshenry/cobrakbase.git && \
    cd cobrakbase && git checkout 68444e46fe3b68482da80798642461af2605e349
RUN mkdir test

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]