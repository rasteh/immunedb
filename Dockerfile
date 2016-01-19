FROM ubuntu:14.04
RUN apt-get update
RUN sudo apt-get install -y python-numpy python-scipy python-setuptools supervisor mysql-server
COPY setup.py /app/
COPY sldb/ /app/sldb
COPY lib/ /app/lib
COPY bin/ /app/bin
WORKDIR /app
RUN ls
RUN python setup.py install
COPY docker/supervisord.conf /etc/supervisor
COPY docker/my.cnf /etc/mysql
CMD ["/usr/bin/supervisord", "-c", "/etc/supervisor/supervisord.conf"]
