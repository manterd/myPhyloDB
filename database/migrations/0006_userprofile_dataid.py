# -*- coding: utf-8 -*-
# Generated by Django 1.9.2 on 2016-09-24 02:09
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('database', '0005_field_changes'),
    ]

    operations = [
        migrations.AddField(
            model_name='userprofile',
            name='dataID',
            field=models.CharField(blank=True, max_length=200),
        ),
    ]
