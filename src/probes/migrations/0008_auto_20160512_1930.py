# -*- coding: utf-8 -*-
# Generated by Django 1.9.1 on 2016-05-12 11:30
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('probes', '0007_sample_offset'),
    ]

    operations = [
        migrations.AlterField(
            model_name='cellline',
            name='primary_hist',
            field=models.CharField(max_length=100),
        ),
    ]
