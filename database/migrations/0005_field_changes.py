# -*- coding: utf-8 -*-
# Generated by Django 1.9.2 on 2016-08-05 04:51
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('database', '0004_userprofile'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='reference',
            options={},
        ),
        migrations.AlterField(
            model_name='project',
            name='end_date',
            field=models.CharField(blank=True, max_length=25),
        ),
        migrations.AlterField(
            model_name='project',
            name='pi_affiliation',
            field=models.CharField(blank=True, max_length=100),
        ),
        migrations.AlterField(
            model_name='project',
            name='pi_first',
            field=models.CharField(blank=True, max_length=50),
        ),
        migrations.AlterField(
            model_name='project',
            name='pi_last',
            field=models.CharField(blank=True, max_length=75),
        ),
        migrations.AlterField(
            model_name='project',
            name='project_name',
            field=models.CharField(blank=True, max_length=100),
        ),
        migrations.AlterField(
            model_name='project',
            name='projectid',
            field=models.CharField(max_length=50, primary_key=True, serialize=False),
        ),
        migrations.AlterField(
            model_name='project',
            name='start_date',
            field=models.CharField(blank=True, max_length=25),
        ),
        migrations.AlterField(
            model_name='reference',
            name='refid',
            field=models.CharField(max_length=50, primary_key=True, serialize=False),
        ),
    ]
