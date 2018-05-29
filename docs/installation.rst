.. highlight:: shell

============
Installation
============

SMETANA currently supports Python 2.7 and Python 3.6, which are available for all major operating systems. We recommend the `Anaconda python
distribution <https://www.continuum.io/downloads>`_.

It can be easily installed using the **pip** package manager:

.. code-block:: console

    $ pip install smetana

This will automatically install other dependencies as well:

- framed_
- pandas_

.. _framed: https://github.com/cdanielmachado/framed
.. _pandas: https://pandas.pydata.org/

Additionally, you must install one of the supported solvers:

- CPLEX_
- Gurobi_

.. _CPLEX: https://www.ibm.com/analytics/data-science/prescriptive-analytics/cplex-optimizer
.. _Gurobi: https://www.gurobi.com/

From our experience you get faster results using CPLEX.
