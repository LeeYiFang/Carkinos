# -*- coding: utf-8 -*-
# Generated by Django 1.9 on 2016-01-06 15:07
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('probes', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='ProbeID',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('Probe_id', models.CharField(max_length=20)),
                ('Gene_symbol', models.CharField(max_length=20)),
                ('Entrez_id', models.IntegerField()),
                ('Gene_name', models.TextField(blank=True, default='')),
            ],
        ),
        migrations.DeleteModel(
            name='ProbeID_GeneSymbol',
        ),
    ]
