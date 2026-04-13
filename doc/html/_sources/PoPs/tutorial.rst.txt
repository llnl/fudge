.. _tutorial:

Tutorial and Examples
=====================

This tutorial provides a quick start, showing how to create a PoPs database, populate it with some
particles and then access particle properties. Start by creating a new empty database:

>>> from PoPs import database
>>> from PoPs.families import baryon, nuclide
>>> from PoPs.quantities import quantity
>>> pops = database.database(name="PoPs_example", version="1.0")
>>> defaultLabel = 'eval'

The defaultLabel is used to identify particle properties (see discussion about multiple assignments in
the `quantities module <quantities.html>`_).
Next we create a new particle (the neutron), add some properties and then add it to the database:

>>> neutron = baryon.particle('n')
>>> neutron.buildFromRawData(
    mass = (1.00866491574, 'amu'),
    spin = (1/2, 'hbar'),
    parity = (1,''),
    charge = (0,'e'),
    halflife = (881.5,'s'),
    label = defaultLabel
)
>>> pops.add(neutron)

The 'add' method determines where the new particle should be stored in the database (in this case it belongs
in the 'baryons' section) and stores it appropriately.

Uncertainties can also be added to a quantity:

>>> from xData.uncertainty.physicalQuantity import uncertainty, standard
>>> neutron.mass[defaultLabel].uncertainty = uncertainty.uncertainty(
    standard.standard( uncertainty.double(4.9e-10) ) )

.. note:: Helper methods may be introduced in the future to facilitate adding uncertainty.

The uncertainty is assumed to have the same units as its associated value.

A PoPs particle can print itself to XML:

>>> print( neutron.toXML() )
<baryon id="n">
  <mass>
    <double label="eval" value="1.00866491574" unit="amu">
      <uncertainty>
        <standard>
          <double value="1.5e-08"/></standard></uncertainty></double></mass>
  <spin>
    <fraction label="eval" value="0" unit="hbar"/></spin>
  <parity>
    <integer label="eval" value="1"/></parity>
  <charge>
    <integer label="eval" value="0" unit="e"/></charge>
  <halflife>
    <double label="eval" value="881.5" unit="s"/></halflife></baryon>

Next we add an excited-state nuclide to the database:

>>> Am242_e2 = nuclide.particle( 'Am242_e2' )

When a nuclide (i.e. an atom containing a particular nucleus) is created, the corresponding nucleus
is also automatically created. Properties can be added to both the nuclide and the nucleus, where
the nuclide properties apply to the electrons as well as the nucleus.

>>> Am242_e2.buildFromRawData( mass=(242.0593557, 'amu') )  # atomic mass including electrons
>>> Am242_e2.nucleus.buildFromRawData( spin=(5, 'hbar') )   # nuclear spin, not including electrons

The buildFromRawData method does not yet support adding excited state energies, but they can be added directly:

>>> from PoPs.quantities import nuclearEnergyLevel
>>> Am242_e2.nucleus.energy.add( nuclearEnergyLevel.double(4.86e+4, 'eV') )

Now add this new particle to the database. Note that this will also add a chemicalElement particle group
(Americium) and an isotope particle group (Am242), inserting the nuclide inside the isotope.

>>> pops.add( Am242_e2 )

You may recognize Am242_e2 as a metastable state. We can add an alias identifying it as metastable:

>>> from PoPs.alias import MetaStable
>>> pops.aliases.add( MetaStable( id='Am242_m1', pid=Am242_e2.id, metaStableIndex=1) )

This gives us two different ways to access the metastable state:

>>> tmp1 = pops['Am242_e2']
>>> tmp2 = pops['Am242_m1']
>>> tmp2.pid == tmp1.id
True

A PoPs database can be written out to an XML file:

>>> pops.saveToFile("pops.xml")

To load the data back into memory,

>>> popsCopy = database.readFile("pops.xml")
