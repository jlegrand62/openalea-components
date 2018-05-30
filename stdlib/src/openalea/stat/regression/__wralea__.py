# -*- python -*-
#
#       OpenAlea.StdLib
#
#       Copyright 2006 - 2008 INRIA - CIRAD - INRA  
#
#       File author(s): CHAUBERT Florence <florence.chaubert@cirad.fr>
#                       Da SILVA David <david.da_silva@cirad.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
# 
#       OpenAlea WebSite : http://openalea.gforge.inria.fr
#


from openalea.core import Factory as Fa
from openalea.core import IBool, IDict, IEnumStr, IFloat, ISequence

__name__ = "openalea.stat.regression"
__alias__ = ["stat.regression"]

__version__ = '0.0.1'
__license__ = 'CECILL-C'
__authors__ = 'OpenAlea consortium'
__institutes__ = 'INRIA/CIRAD/UM2'
__description__ = 'Regressions from Rpy and Scipy.'
__url__ = 'http://rpy.sourceforge.net and http://www.scipy.org/'

__editable__ = 'False'

__all__ = ['glm', 'linearregression', 'multiplelinearregression',
           'linearregressiontoplot', 'linearregressionscipy', ]

glm = Fa(uid="e7f010b84e7811e6bff6d4bed973e64a",
         name="glm (rpy)",
         description="Compute the generalized linear regression",
         category="openalea.stat.regression.regression",
         nodemodule="regression",
         nodeclass="Glm",
         inputs=(dict(name="X", interface=ISequence, showwidget=True),
                 dict(name="Y", interface=ISequence, showwidget=True),
                 dict(name="Family", interface=IEnumStr(
                     ['binomial', 'Gamma', 'gaussian', 'poisson']),
                      showwidget=True),
                 ),
         outputs=(dict(name="Glm", interface=IDict),
                  ),
         )

linearregression = Fa(uid="e7f010b94e7811e6bff6d4bed973e64a",
                      name="linear regression (rpy)",
                      description="Compute the linear regression",
                      category="regression",
                      nodemodule="openalea.stat.regression.regression",
                      nodeclass="LinearRegression",
                      inputs=(
                          dict(name="X", interface=ISequence, showwidget=True),
                          dict(name="Y", interface=ISequence, showwidget=True),
                          dict(name="alpha", interface=IFloat, value=5.),
                          dict(name="origin", interface=IBool, value=False),
                      ),
                      outputs=(dict(name="reg", interface=IDict),
                               ),
                      )

multiplelinearregression = Fa(uid="e7f010ba4e7811e6bff6d4bed973e64a",
                              name="multiple linear regression (rpy)",
                              description="Compute a multiple linear regression",
                              category="regression",
                              nodemodule="openalea.stat.regression.regression",
                              nodeclass="multiReg",
                              inputs=(dict(name="X", interface=ISequence,
                                           showwidget=True),
                                      dict(name="Y", interface=ISequence,
                                           showwidget=True),
                                      dict(name="colList", interface=ISequence,
                                           showwidget=True),
                                      dict(name="alpha", interface=IFloat,
                                           value=5.),
                                      ),
                              outputs=(dict(name="reg", interface=IDict),
                                       ),
                              )

linearregressiontoplot = Fa(uid="e7f010bb4e7811e6bff6d4bed973e64a",
                            name="linear regression to plot (plotools)",
                            description="Generate plotable object from linear regression",
                            category="regression",
                            nodemodule="openalea.stat.regression.regression",
                            nodeclass="LR2Plot",
                            inputs=(dict(name='reg', interface=IDict),
                                    ),
                            outputs=(dict(name='plotObjList0', ),
                                     dict(name='plotObjList1', ),
                                     ),
                            )

linearregressionscipy = Fa(uid="e7f010bc4e7811e6bff6d4bed973e64a",
                           name="linear regression (scipy)",
                           description="Compute the linear regression with scipy",
                           category="regression",
                           nodemodule="openalea.stat.regression.regression",
                           nodeclass="linearregress",
                           inputs=(
                               dict(name="X", interface=ISequence,
                                    showwidget=True),
                               dict(name="Y", interface=ISequence,
                                    showwidget=True),
                           ),
                           outputs=(dict(name="linearregress", interface=IDict),
                                    ),
                           )
