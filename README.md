# CellExpress

The website is freely available at http://cellexpress.cgm.ntu.edu.tw/  
CellExpress is a microarray-based cell line and clinical sample gene expression analysis online system.
The system was written in Python with Django framework, and contains data from Affymetrix U133A and Affy-metrix U133Plus2.0 platforms. 
CellExpress is under the Carkinos project in Bioinformatics and Biostatistics Core, NTU Center of Genomic Medicine.


## Getting Started

### Requirements

- Git
- Python 3.4+
- Django 1.9
- R

#### Binary datasets
All microarrays obtained from the same platform were normalized using the quantile algorithm.The control probes were removed.
Carkinos Project ships several processed datasets including:

- CCLE (GSE36133)
- NCI60 (GSE32474)
- Sanger (GSE68950)
- expO (GSE2109)
- Roth (GSE7307)

The normalized expression data are stored in Numpy binary format, which are generally too large to fit in the GitHub git repository. Make sure you have them first.


### Set up a Python Virtual Environment

#### Built-in `venv`

Create the virtual environment at `<repo>/venv`:

    python3 -m venv venv

To activate it:

    . venv/bin/activate

Or on Windows, use

    . venv\Scripts\activate.bat


### Install Dependencies

Use pip:

    pip install -r requirements.txt


### Set up Local Environment Variables and Database

Settings are stored in environment variables via [django-environ]. The
quickiest way to start is to copy `local.sample.env` into `local.env`:

    cp src/Carkinos/settings/local.sample.env src/Carkinos/settings/local.env

Then edit the `SECRET_KEY` line in `local.env`, replacing `{{ secret_key }}` into any [Django Secret Key] value. An example:

    SECRET_KEY=7h$ab<Lgg!<}|_i\WNxBV^;JUEBC;f~/.-9`9w)G^kqyc:zbYo_s%j]cO5{R!Y=|

After that, follow the instruction to run database migration.


### Go Develop

Change to the `src` directory:

    cd src

Run database migration:

    python manage.py migrate

Run the development server

    python manage.py runserver



## Build Binary Datasets

**NOTE: Binary datasets match the records in SQLite database, make sure they are in sync.**

Please contact s89066[at]gmail[dot]com to let me now if you need to get the binary data in CellExpress.  
For large files, we will provide FTP for you:)

[django-environ]: http://django-environ.readthedocs.org/en/latest/
[Django Secret Key]: http://www.miniwebtool.com/django-secret-key-generator/
