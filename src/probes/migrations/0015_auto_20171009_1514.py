# -*- coding: utf-8 -*-
# Generated by Django 1.9.5 on 2017-10-09 07:14
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('probes', '0014_delete_stand_gene'),
    ]

    operations = [
        migrations.AddField(
            model_name='cellline',
            name='disease',
            field=models.CharField(default='N/A', max_length=100),
        ),
        migrations.AddField(
            model_name='cellline',
            name='organ',
            field=models.CharField(default='N/A', max_length=50),
        ),
    ]
