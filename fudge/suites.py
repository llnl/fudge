# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>


import os
import string
import abc

from LUPY import ancestry as ancestryModule
from fudge import GNDS_formatVersion as GNDS_formatVersionModule

from xData import enums as xDataEnumsModule
from xData import link as linkModule
from PoPs.families import particle as particleModule

class Suite(ancestryModule.AncestryIO_bare, abc.ABC):
    """
    Base class for a class member that is list like. For example, the lists inside the class 
    reactionSuite ('reactions', 'sums', 'productions' and 'fissionComponents').
    """

    def __init__( self, allowedClasses, allow_href = False ) :

        ancestryModule.AncestryIO_bare.__init__(self)
        self.__allowedClasses = [ cls for cls in allowedClasses ]

        self.__items = []
        self.__allow_href = allow_href
        self.__href = None
        self.__hrefInstance = None

    def __contains__( self, label ) :

        hrefInstance = self.hrefInstance( )
        if( hrefInstance is not None ) : return( label in hrefInstance )

        for item in self.__items :
            if( item.label == label ) : return( True )
        return( False )

    def __delitem__( self, key ) :

        if( self.hrefInstance( ) is not None ) : raise AttributeError( 'href suite does not support modification of href instance.' )

        if( isinstance( key, int ) ) :
            del self.__items[key]
            return
        status = self.remove( key )
        if( not( status ) ) : KeyError( 'key not in suite' )

    def __getitem__( self, label ) :

        hrefInstance = self.hrefInstance( )
        if( hrefInstance is not None ) : return( hrefInstance[label] )

        if isinstance(label, slice): return [self.__items[index] for index in range(*label.indices(len(self)))]
        if isinstance(label, int): return self.__items[label]

        if( not( isinstance( label, str ) ) ) : raise TypeError( "label must be a string" )
        for item in self.__items :
            if( item.label == label ) : return( item )

        # requested style not found, but what about styles it derives from?
        if hasattr(self.rootAncestor, 'styles'):
            requestedStyle = self.rootAncestor.styles[ label ]
            parentStyle = requestedStyle.derivedFromStyle

            while parentStyle is not None:
                if parentStyle.label in self:
                    return self[parentStyle.label]
                parentStyle = parentStyle.derivedFromStyle

        raise KeyError( "item with label '%s' not found in suite '%s'" % ( label, self.moniker ) )

    def __iter__( self ) :

        hrefInstance = self.hrefInstance( )
        if( hrefInstance is not None ) :
            hrefInstance.__iter__( )
        else :
            n1 = len( self )
            for i1 in range( n1 ) : yield self.__items[i1]

    def __len__( self ) :

        hrefInstance = self.hrefInstance( )
        if( hrefInstance is not None ) : return( len( hrefInstance ) )

        return( len( self.__items ) )

    @property
    def allowedClasses( self ) :

        return( self.__allowedClasses )

    @property
    def href( self ) :

        return( self.__href )

    @property
    def allow_href( self ) :

        return( self.__allow_href )

    def add(self, newItem, setAncestor = True, addLazyParsingHelper=False):
        """
        Adds newItem to the suite. If another item in the suite has the same label as newItem, a KeyError is raised.
        If newItem is not an allowed class, a TypeError is raised.

        :param newItem:
        :return:
        """

        if( self.hrefInstance( ) is not None ) : raise AttributeError( 'href suite does not support modification of href instance.' )

        if( not( isinstance( newItem.label, str ) ) ) : raise IndexError( '''Item's label must be a string: %s''' % type( newItem.label ) )

        allowedClasses = [ allowedClass for allowedClass in self.__allowedClasses ]
        if hasattr(self, 'specialAllowedClasses'):                      # Currently, this happends in incompleteReactions.
            for allowedClass in getattr(self, 'specialAllowedClasses'): allowedClasses.append(allowedClass)

        if not addLazyParsingHelper:
            found = False
            for cls in allowedClasses :
                if( isinstance( newItem, cls ) ) :
                    found = True
                    break
            if not found: raise TypeError( 'Invalid class "%s" for suite "%s"' % ( newItem.__class__, self.moniker ) )

        for i1, item in enumerate( self.__items ) :
            if( item.label == newItem.label ) : raise KeyError( 'item with label = "%s" already present in suite: %s.' % (item.label, self.toXLink()))

        if( setAncestor ) : newItem.setAncestor(self)
        self.__items.append( newItem )

    def amendForPatch( self, fromLabel, toLabel ) :

        if( self.hrefInstance( ) is None ) :
            for item in self: item.amendForPatch(fromLabel, toLabel)

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        if( self.hrefInstance( ) is None ) :
            for index, item in enumerate(self): item.convertUnits(unitMap)

    def checkAncestry( self, verbose = 0, level = 0 ) :

        if( self.hrefInstance( ) is None ) :
            for item in self: item.checkAncestry(verbose = verbose, level = level)

    def diff( self, other, diffResults ) :

        if( self.hrefInstance( ) is not None ) : return

        labels1 = [ item.label for item in self ]
        labels2 = [ item.label for item in other ]
        labels = [ label for label in labels1 ]
        for label in labels2 :
            if( label not in labels ) : labels.append( label )

        for label in labels :
            if( label not in labels1 ) :
                diffResults.append( '%s missing - 1' % other.moniker, '', '', other[label].toXLink( ) )
            elif( label not in labels2 ) :
                diffResults.append( '%s missing - 2' % self.moniker, '', self[label].toXLink( ), '' )
            else :
                if( type( self ) != type( other ) ) :
                    print( '    Cannot diff type "%s" with type "%s".' % ( type( self ), type( other ) ) )
                else :
                    if( hasattr( self[label], 'diff' ) ) :
                        self[label].diff( other[label], diffResults )
                    else :
                        print( '    Method "diff" missing for object %s' % type( self[label] ) )

    def findInstancesOfClassInChildren( self, cls, level = 9999 ) :

        foundInstances = []
        level -= 1
        if( level < 0 ) : return( foundInstances )
        for item in self :
            if( isinstance( item, cls ) ) : foundInstances.append( item )
            foundInstances += item.findInstancesOfClassInChildren( cls, level )

        return( foundInstances )

    def findLinks( self, links ) :

        if( self.hrefInstance( ) is not None ) : return

        for item in self :
            if( hasattr( item, 'findLinks' ) ) :
                getattr( item, 'findLinks' )( links )
            else :
                if( isinstance( item, ( linkModule.Link, linkModule.Link2 ) ) ) : links.append( [ item, item.link, item.path ] )

    def fixDomains(self, labels, energyMin, energyMax):
        """
        Class **fixDomains** for each entry of *self* if it has a **fixDomains** method.
        """

        numberOfFixes = 0
        for entry in self:
            if entry.label in labels:
                if hasattr(entry, 'fixDomains'):
                    numberOfFixes += entry.fixDomains(energyMin, energyMax, xDataEnumsModule.FixDomain.both)
                else:
                    print('WARNING: type "%s" does not have a "fixDomains" method.' % type(entry))

        return numberOfFixes

    def hrefInstance( self ) :

        if( self.__hrefInstance is None ) :
            if( self.__allow_href ) :
                if( self.__href is not None ) :
                    self.__hrefInstance = self.rootAncestor.followXPath( self.__href )

        return( self.__hrefInstance )

    def set_href( self, href ) :

        if( not( self.__allow_href ) ) : raise Exception( 'Suite does not allow href.' )
        if( len( self.__items ) != 0 ) : raise Exception( 'Suite contains items and cannot also have an href.' )
        if( not( isinstance( href, str ) ) ) : raise ValueError( 'href instance must be a string.' )

        self.__href = href

    def pop( self, label, *args ) :
        """
        Remove item by label, and return it. If the label is not found, raise KeyError
        :param label:
        :param args:
        :return:
        """

        if( self.hrefInstance( ) is not None ) : raise AttributeError( 'href suite does not support modification of href instance.' )

        if( len( args ) > 2 ) : raise Exception( 'Only one default value is allowed: got %s' % len( args ) )
        for i1, item in enumerate( self.__items ) :
            if( item.label == label ) : return( self.__items.pop( i1 ) )
        if( len( args ) == 1 ) : return( args[0] )
        raise KeyError( label )

    def cullStyles( self, styleList ) :

        if( self.hrefInstance( ) is None ) :
            for item in self: item.cullStyles(styleList)

    def labels( self ) :
        """
        Returns the list of all the labels in self.
        """

        hrefInstance = self.hrefInstance( )
        if( hrefInstance is not None ) : return( hrefInstance.labels( ) )

        return( [ item.label for item in self.__items ] )

    def printAllowedChildren(self):
        """
        This method prints all the allowed classes that can be a added as a child of *self*.
        """

        for cls in self.allowedClasses:
            try:
                location = str(cls).split("'")[1]
            except:
                location = str(cls)
            print('    %s' % location)

    def remove( self, label ) :
        """
        Remove item by label. Returns True if label was present, otherwise returns False
        :param label: str
        :return: bool
        """

        if( self.hrefInstance( ) is not None ) : raise AttributeError( 'href suite does not support modification of href instance.' )

        for i1, item in enumerate( self.__items ) :
            if( item.label == label ) :
                del self.__items[i1]
                return( True )
        return( False )

    def removeStyles( self, styleLabels ) :
        """
        Removes all forms whose label matches one of the labels in removeStyles.
        """

        if( self.hrefInstance( ) is not None ) : raise AttributeError( 'href suite does not support modification of href instance.' )

        for item in self : item.removeStyles( styleLabels )

    def replace( self, newItem ) :
        """
        Replace an existing item with newItem. If no item in the suite has the same label as newItem, a KeyError is raised.
        If newItem is not an allowed class, a TypeError is raised.

        :param newItem:
        :return:
        """

        if( self.hrefInstance( ) is not None ) : raise AttributeError( 'href suite does not support modification of href instance.' )

        if( not( isinstance( newItem.label, str ) ) ) : raise ValueError( '''Item's label must be a string''' )
        found = False
        for cls in self.__allowedClasses :
            if( isinstance( newItem, cls ) ) :
                found = True
                break
        if( not( found ) ) : raise TypeError( 'Invalid class "%s" for suite "%s"' % ( newItem.__class__, self.moniker ) )

        for i1, item in enumerate( self.__items ) :
            if( item.label == newItem.label ) :
                del self.__items[i1]
                self.__items.insert( i1, newItem )
                newItem.setAncestor(self)
                return

        raise KeyError( 'item with label = "%s" not present in suite' % newItem.label )

    def clear( self ):
        """
        Remove all members of the suite
        :return:
        """

        if( self.hrefInstance( ) is not None ) : raise AttributeError( 'href suite does not support modification of href instance.' )

        self.__items = []

    def toXML_strList(self, indent='', **kwargs):

        if self.__href is not None:
            return ['%s<%s href="%s"/>' % (indent, self.moniker, self.__href)]

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        moniker = self.monikerByFormat.get(kwargs.get('formatVersion'), self.moniker)
        if len(self) == 0:
            showEmptySuite = kwargs.get('showEmpty', False) or kwargs.get('showEmptySuite', False)
            if showEmptySuite:
                return ['%s%-24s <!-- suite -->' % (indent, ('<%s/>' % moniker))]
            return []

        xmlString = ['%s<%s>' % (indent, moniker)]
        for item in self:
            xmlString += item.toXML_strList(indent2, **kwargs)
        xmlString[-1] += '</%s>' % moniker

        return xmlString

    def uniqueLabel( self, item ) :
        """
        If item's label is the same as another item's label in self, construct a new unique label
        based on item's label appended with '__' and one or more lower case letters (i.e., 'a' to 'z').
        """

        def incrementSuffix( suffix ) :

            if( len( suffix ) == 0 ) : return( 'a' )
            index = string.ascii_lowercase.find( suffix[-1] ) + 1
            if( index != 26 ) : return( suffix[:-1] + string.ascii_lowercase[index] )
            return( incrementSuffix( suffix[:-1] ) + 'a' )

        if( item.label in self ) :
            label = item.label
            label__ = label + '__'
            n1 = len( label__ )
            l1 = 0
            suffixes = []
            for _item in self.__items :          # Find list of longest labels that start with label__.
                if( _item.label[:n1] == label__ ) :
                    suffix = _item.label[n1:]
                    if( not( suffix.islower( ) ) ) : continue       # Ignore non-standard labels.
                    l2 = len( suffix )
                    if( l2 < l1 ) : continue
                    if( l2 > l1 ) :
                        l1 = l2
                        suffixes = []
                    suffixes.append( suffix )
            if( len( suffixes ) == 0 ) :
                suffix = 'a'
            else :
                suffix = incrementSuffix( sorted( suffixes )[-1] )
            item.label = label__ + suffix
        return( item )

    def parseNode(self, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)

        allowedClasses = {}
        for allowedClass in self.allowedClasses: allowedClasses[allowedClass.moniker] = allowedClass

        nodesNotParsed = []
        for child in node:
            cls = allowedClasses.get(child.tag)
            if cls is None:
                cls = allowedClasses.get(self.legacyMemberNameMapping[child.tag])
                if cls is None:
                    nodesNotParsed.append(child.tag)
                    continue
            self.add(cls.parseNodeUsingClass(child, xPath, linkData, **kwargs))

        if len(nodesNotParsed) > 0: raise Exception("Encountered unexpected child nodes '%s'!" % ', '.join(nodesNotParsed))

        xPath.pop()

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        This method must be overrwritten by the derived class.
        """

        xPath.append(node.tag)
        instance = cls( )
        xPath.pop()

        instance.parseNode(node, xPath, linkData, **kwargs)

        return instance

class ExclusiveSuite(Suite):

    def fixDomains(self, labels, energyMin, energyMax):
        """
        Calls the **fixDomains** for each entry in self.
        """

        numberOfFixes = 0
        for entry in self: numberOfFixes += entry.fixDomains(labels, energyMin, energyMax)

        return numberOfFixes

class ExternalFiles(ExclusiveSuite):

    moniker = 'externalFiles'

    def __init__( self ) :

        from fudge import externalFile as externalFileModule

        ExclusiveSuite.__init__( self, [ externalFileModule.ExternalFile ] )

class Reactions(ExclusiveSuite):

    moniker = 'reactions'

    def __init__( self ) :

        from fudge.reactions import reaction as reactionModule

        ExclusiveSuite.__init__( self, [ reactionModule.Reaction ] )

    def asSortedList( self, PoPs ) :
        """Returns a list of the reactions sorted. This method uses a reaction's reactionProducts method to determine sorting order."""

        class Particle( particleModule.Particle ) :

            moniker = 'dummy'
            familyOrder = 999

            def __lt__( self, other ) :

                if( self.familyOrderLessThan( other ) ) : return( True )        # Also checks isinstance of other.
                if( self.familyOrder != other.familyOrder ) : return( False )
                return( self.id < other.id )

        listToSort = []
        for index, reaction in enumerate( self ) :
            particles, intermediates, processes, _reaction = reaction.reactionProducts( ).asSortedList( PoPs, excludePhotons = True )
            items = [ [ PoPs[pid], value ] for pid, value in particles + intermediates ]
            for item in items :
                if( item[1] == 0 ) : item[1] = 999                          # Make a large number so its at the end of sorting if all else equal.
            items.append( [ Particle( key[0] ) for key in processes ] )
            items.append( Particle( str( index ) ) )
            items[-1].reaction = reaction
            listToSort.append( items )
        listToSort.sort( )

        return( [ item[-1].reaction for item in listToSort ] )

    def diff( self, other, diffResults ) :

        class TwoPoPs :

            def __init__( self, pops1, pops2 ) :

                self.pops1 = pops1
                self.pops2 = pops2

            def __getitem__( self, key ) :

                try :
                    return( self.pops1[key] )
                except :
                    return( self.pops2[key] )

        sortedReactions1 = self.asSortedList( self.ancestor.PoPs )
        sortedReactions2 = other.asSortedList( other.ancestor.PoPs )
        twoPoPs = TwoPoPs( self.ancestor.PoPs, other.ancestor.PoPs )

        n1 = len( sortedReactions1 )
        index1 = 0
        n2 = len( sortedReactions2 )
        index2 = 0
        reactionPairs = []
        while( ( index1 < n1 ) and ( index2 < n2 ) ) :
            reaction1 = sortedReactions1[index1]
            reaction2 = sortedReactions2[index2]
            if( reaction1.label == reaction2.label ) :
                compare = '='
            else :
                if( reaction1.reactionProducts( ).cullPhotons( ) == reaction2.reactionProducts( ).cullPhotons( ) ) :
                    compare = '='
                else :
                    dummy = Reactions( )
                    dummy.add( reaction1, setAncestor = False )
                    dummy.add( reaction2, setAncestor = False )
                    dummy = dummy.asSortedList( twoPoPs )
                    compare = '<'
                    if( dummy[0] is reaction2 ) : compare = '>'

            if( compare == '=' ) :
                reactionPairs.append( [ reaction1, reaction2 ] )
                index1 += 1
                index2 += 1
            elif( compare == '<' ) :
                reactionPairs.append( [ reaction1, None ] )
                index1 += 1
            else :
                reactionPairs.append( [ None, reaction2 ] )
                index2 += 1

        for reaction1, reaction2 in reactionPairs :
            if( reaction1 is None ) :
                diffResults.append( '%s missing - 1' % other.moniker, '', '', other[reaction2.label].toXLink( ) )
            elif( reaction2 is None ) :
                diffResults.append( '%s missing - 2' % self.moniker, '', self[reaction1.label].toXLink( ), '' )
            else :
                reaction1.diff( reaction2, diffResults )

    def reactionProducts( self, PoPs = None, excludePhotons = True ) :

        if( PoPs is None ) : PoPs = self.findAttributeInAncestry( 'PoPs' )

        sortedReactions = self.asSortedList( PoPs )
        return( [ reaction.reactionProducts( ) for reaction in sortedReactions ] )

class OrphanProducts(ExclusiveSuite):

    moniker = 'orphanProducts'

    def __init__( self ) :

        from fudge.reactions import orphanProduct as orphanProductModule

        ExclusiveSuite.__init__( self, [ orphanProductModule.OrphanProduct ] )

class CrossSectionSums(ExclusiveSuite):

    moniker = 'crossSectionSums'
    monikerByFormat = { GNDS_formatVersionModule.version_1_10: 'crossSections', GNDS_formatVersionModule.version_2_0_LLNL_3 : 'crossSections' }

    def __init__( self ) :

        from fudge import sums as sumsModule

        ExclusiveSuite.__init__( self, [sumsModule.CrossSectionSum] )

class MultiplicitySums(ExclusiveSuite):

    moniker = 'multiplicitySums'
    monikerByFormat = { GNDS_formatVersionModule.version_1_10: 'multiplicities', GNDS_formatVersionModule.version_2_0_LLNL_3 : 'multiplicities' }

    def __init__( self ) :

        from fudge import sums as sumsModule

        ExclusiveSuite.__init__( self, [sumsModule.MultiplicitySum] )

class Productions(ExclusiveSuite):

    moniker = 'productions'

    def __init__( self ) :

        from fudge.reactions import production as productionModule

        ExclusiveSuite.__init__( self, [productionModule.Production] )

class IncompleteReactions(ExclusiveSuite):

    moniker = 'incompleteReactions'

    def __init__( self ):

        from .reactions import reaction as reactionModule
        from .reactions import incompleteReaction as incompleteReactionModule

        ExclusiveSuite.__init__(self, [incompleteReactionModule.IncompleteReaction])
        self.specialAllowedClasses = (reactionModule.Reaction, )                            # Allow reaction as it is needed at times but do not remember why.

class FissionComponents(ExclusiveSuite):

    moniker = 'fissionComponents'

    def __init__( self ) :

        from fudge.reactions import fissionComponent as fissionComponentModule

        ExclusiveSuite.__init__( self, [fissionComponentModule.FissionComponent] )

class ApplicationData(Suite):

    moniker = 'applicationData'

    def __init__( self ):

        from fudge import institution as institutionModule

        Suite.__init__( self, [ institutionModule.Institution, institutionModule.UnknownInstitution ] )
