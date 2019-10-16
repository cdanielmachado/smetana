#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from smetana.interface import main
import pandas as pd


class TestGlobal(unittest.TestCase):
    def test_global(self):
        main(["tests/data/ec_*_ko.xml"], mode="global", output="tests/output/test", media="M9,LB",
             mediadb="tests/data/media_db.tsv", ext_comp_id="e", exclude="tests/data/inorganic.txt")
        df = pd.read_csv("tests/output/test_global.tsv")
        self.assertEqual(df.shape[0], 2)


class TestDetailed(unittest.TestCase):

    def test_detailed(self):
        main(["tests/data/ec_*.xml"], mode="detailed", output="tests/output/test", media="M9,LB",
             mediadb="tests/data/media_db.tsv", ext_comp_id="e", exclude="tests/data/inorganic.txt")
        df = pd.read_csv("tests/output/test_detailed.tsv")
        self.assertGreater(df.shape[0], 5)
        self.assertLess(df.shape[0], 15)
