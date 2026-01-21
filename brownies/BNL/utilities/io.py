import os.path


def read_evaluation(endfFile, reconstructResonances=True):
    """
    Read in the evaluation
    :param endfFile: the file name.  This routine will attempt to determine the file type itself.
    :param reconstructResonances: FIXME: add documentation
    :return: the evaluation as a dict of results like {'reactionSuite':the reaction suite, 'covarianceSuite':
             the covariances}
    """
    if not os.path.exists(endfFile):
        raise IOError("File named " + endfFile + " doesn't exist!")
    with open(endfFile) as evaluationfile:
        firstline = evaluationfile.readline()
    if firstline.startswith("<?xml") or firstline.startswith("<reactionSuite "):
        from fudge import reactionSuite
        return {'reactionSuite': reactionSuite.ReactionSuite.readXML_file(endfFile), 'covarianceSuite': None, 'info': {}, 'errors': []}
    else:
        import brownies.legacy.toENDF6.toENDF6     # this import adds 'toENDF6' methods to many GNDS classes
        from brownies.legacy.converting import endfFileToGNDS
        return endfFileToGNDS.endfFileToGNDS(endfFile, toStdOut=False, skipBadData=True, continuumSpectraFix=True,
                                             reconstructResonances=reconstructResonances)
