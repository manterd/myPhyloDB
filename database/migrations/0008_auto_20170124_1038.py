# -*- coding: utf-8 -*-
# Generated by Django 1.9.5 on 2017-01-24 17:38
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('database', '0007_auto_20170120_1533'),
    ]

    operations = [
        migrations.AlterField(
            model_name='userprofile',
            name='dataID',
            field=models.CharField(blank=True, default=b'TEMPID', max_length=200),
        ),
    ]