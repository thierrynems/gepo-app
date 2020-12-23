FROM python:3.6
#permet à python d'afficher tout ce aui est imprimé dans l'application plutot que de le mettre en memoire tampon.
ENV PYTHONUNBUFFERED 1

#creation du repertoire de travail / WORKDIR
RUN mkdir /app
WORKDIR /app
#ajout du fichier requirements.txt dans /app
ADD requirements.txt /app

#installation requirements
RUN pip install --upgrade pip
RUN pip install -r /app/requirements.txt

