from database.models import Reference

m = Reference.objects.get(projectid='25a0d68af01844949daa2f474be65528')
m.alignDB = 'silva.seed_v119.align'
m.taxonomyDB = 'gg_13_5_99.pds.tax'
m.templateDB = 'gg_13_5_99.fasta'
m.save()