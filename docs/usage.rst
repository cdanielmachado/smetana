=====
Usage
=====


Single community
________________


To run SMETANA on a singe community just give it a list of single-species models in SBML format:

.. code-block:: console

    $ smetana species1.xml species2.xml ...

    $ smetana *.xml

This will generate two files called ``global.tsv`` and ``detailed.tsv`` with global and detailed properties on the
given microbial community (see Algorithms section for details).


Multiple communities
____________________

You can run SMETANA on multiple communities in a single line of code. For that you need to create a table describing the
composition of each community. This table is expected to be in long format, tab-separated, with two columns:
**community id**, **organism id**. The latter should match the SBML filename (without extension). Example:

+----------+---------+
|community1|organism1|
+----------+---------+
|community1|organism2|
+----------+---------+
|community2|organism1|
+----------+---------+
|community2|organism3|
+----------+---------+

After creating this file, invoke SMETANA as follows:

.. code-block:: console

    $ smetana *.xml -c communities.tsv

Medium composition
__________________

By default, SMETANA will simulate every community using a complete medium composition. To test your communities on
different media, give it the list of media and a file containing the media compositions. You can specify multiple media
in a single command and SMETANA will run consecutively for each medium.

.. code-block:: console

    $ smetana *.xml -m M9,LB --mediadb library.tsv

For an example on how to create your own library file please check this example_.

.. _example: https://github.com/cdanielmachado/carveme/blob/master/carveme/data/input/media_db.tsv


Advanced
________

SMETANA offers other options, such as:

- Selecting different SBML *flavors* (``--flavor``)
- Changing the name/directory of your output files (``-o``, ``--output``)
- Selecting only a subset of scores to run (``-s``, ``--scores``)
- Changing your default solver (``--solver``)
- Specifying the identifier of the extracellular compartment in the models (``--ext``).

For more detailed instructions just type:

.. code-block:: console

    $ smetana -h

