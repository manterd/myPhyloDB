# -*- coding: utf-8 -*-
# Generated by Django 1.11.12 on 2018-07-06 15:12
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('database', '0017_otu_99_otuseq'),
    ]

    operations = [
        migrations.AddField(
            model_name='userprofile',
            name='gavePermsTo',
            field=models.TextField(blank=True),
        ),
        migrations.AddField(
            model_name='userprofile',
            name='hasPermsFrom',
            field=models.TextField(blank=True),
        ),
    ]
