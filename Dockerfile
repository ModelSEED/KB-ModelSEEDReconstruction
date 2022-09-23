FROM kbase/sdkbase2:python
MAINTAINER chenry@anl.gov
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# RUN apt-get update

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


RUN rm -rf /miniconda/lib/python3.6/site-packages/numpy
RUN rm -rf /miniconda/lib/python3.6/site-packages/ruamel*
    
# Install forked version of optlang and cobrapy to add
# additional solver support COINOR-CBC,CLP and OSQP.
# Must install cobrakbase first since it installs a newer version
# of cobra not supported on the custom branches for optlang/cobra.
# Running with --ignore-installed will overwrite with the correct
# cobra version.
RUN mkdir deps && cd deps && \
    git clone --branch cobra-model https://github.com/fxe/cobrakbase.git && \
    pip install cobrakbase/ --ignore-installed && \
    git clone https://github.com/ModelSEED/ModelSEEDpy.git && \
    pip install Jinja2 

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]