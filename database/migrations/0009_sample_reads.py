# -*- coding: utf-8 -*-
# Generated by Django 1.9.5 on 2017-02-07 00:57
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('database', '0008_auto_20170124_1038'),
    ]

    operations = [
        migrations.AddField(
            model_name='sample',
            name='reads',
            field=models.IntegerField(default=0, editable=False, null=True),
        ),
    ]
