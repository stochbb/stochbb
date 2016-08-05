import QtQuick 2.1
import QtQuick.Controls 1.4
import QtQuick.Layouts 1.3
//import ToolBarSeparator 1.0
import io.github.stochbb 1.0


ApplicationWindow {
    id: mainwindow
    visible: true
    title: "StochBB"
    minimumWidth: 640
    minimumHeight: 480

    Action {
        id: verifyNet
        text: "Verify"
        iconSource: "qrc:///icons/check_64.png"

    }
    Action {
        id: runNet
        text: "Run"
        iconSource: "qrc:///icons/play_64.png"
    }
    Action {
        id: zoomIn
        text: "Zoom in"
        iconSource: "qrc:///icons/zoom-in_64.png"
    }
    Action {
        id: zoomOut
        text: "Zoom out"
        iconSource: "qrc:///icons/zoom-out_64.png"
    }
    Action {
        id: deleteItem
        text: "Remove"
        iconSource: "qrc:///icons/trash_64.png"
    }
    Action {
        id: showAbout
        text: "About StochBB ..."
    }
    Action {
        id: showManual
        text: "User Manual ..."
    }

    Action {
        id: addStimulus
        text: "Stimulus"
    }
    Action {
        id: addConstant
        text: "Constant"
    }
    Action {
        id: addUniformVariable
        text: "Uniform variable"
    }
    Action {
        id: addNormalVariable
        text: "Normal variable"
    }
    Action {
        id: addGammaVariable
        text: "Gamma variable"
    }
    Action {
        id: addInverseGammaVariable
        text: "inverse Gamma variable"
    }
    Action {
        id: addWeibullVariable
        text: "Weibull variable"
    }
    Action {
        id: addCompoundNormalVariable
        text: "compound Normal variable"
    }
    Action {
        id: addCompoundGammaVariable
        text: "compound Gamma variable"
    }
    Action {
        id: addCompoundInverseGammaVariable
        text: "compound inverse Gamma variable"
    }
    Action {
        id: addCompoundWeibullVariable
        text: "compound Weibull variable"
    }
    Action {
        id: addFixedDelay
        text: "Fixed delay"
    }
    Action {
        id: addGammaProcess
        text: "Gamma process"
    }
    Action {
        id: addInverseGammaProcess
        text: "inverse Gamma process"
    }
    Action {
        id: addWeibullProcess
        text: "Weibull process"
    }
    Action {
        id: addRandomDelay
        text: "Random delay"
    }
    Action {
        id: addCompoundGammaProcess
        text: "compound Gamma process"
    }
    Action {
        id: addCompoundInverseGammaProcess
        text: "compound inverse Gamma process"
    }
    Action {
        id: addCompoundWeibullProcess
        text: "compound Weibull process"
    }
    Action {
        id: addMinimum
        text: "Minimum"
    }
    Action {
        id: addMaximum
        text: "Maximum"
    }
    Action {
        id: addInhibition
        text: "Inhibition"
    }
    Action {
        id: addAffine
        text: "Affine"
    }
    Action {
        id: addMarginalPlot
        text: "Marginal plot"
    }
    Action {
        id: addScatterPlot
        text: "Scatter plot"
    }
    Action {
        id: addKDEPlot
        text: "KDE plot"
    }

    menuBar : MenuBar {
        Menu {
            title: "File"

            MenuItem {
                text: "New"
                shortcut: "Ctrl+N"
                // onTriggered: ...
            }
            MenuItem {
                text: "Open ..."
             shortcut: "Ctrl+O"
            }
            MenuItem {
                text: "Save"
                shortcut: "Ctrl+S"
            }
            MenuItem {
                text: "Save as ..."
                shortcut: "Ctrl+Shift+S"
            }

            MenuSeparator { }

            MenuItem {
                text: "Export network ..."
            }

            MenuSeparator { }

            MenuItem { action: verifyNet }
            MenuItem { action: runNet }

            MenuSeparator { }

            MenuItem {
                text: "Quit"
                shortcut: "Ctrl+Q"
            }
        }

        Menu {
            title: "Edit"

            Menu {
                title: "Add variable"
                iconSource: "qrc:///icons/var_64.png"
                MenuItem { action: addStimulus }
                MenuSeparator {}
                MenuItem { action: addConstant }
                MenuItem { action: addUniformVariable }
                MenuItem { action: addNormalVariable }
                MenuItem { action: addGammaVariable }
                MenuItem { action: addInverseGammaVariable }
                MenuItem { action: addWeibullVariable }
                MenuSeparator {}
                MenuItem { action: addCompoundNormalVariable }
                MenuItem { action: addCompoundGammaVariable }
                MenuItem { action: addCompoundInverseGammaVariable }
                MenuItem { action: addCompoundWeibullVariable }
            }

            Menu {
                title: "Add process"
                iconSource: "qrc:///icons/proc_64.png"
                MenuItem { action: addFixedDelay }
                MenuItem { action: addGammaProcess }
                MenuItem { action: addInverseGammaProcess }
                MenuItem { action: addWeibullProcess }
                MenuSeparator {}
                MenuItem { action: addRandomDelay }
                MenuItem { action: addCompoundGammaProcess }
                MenuItem { action: addCompoundInverseGammaProcess }
                MenuItem { action: addCompoundWeibullProcess }
            }

            Menu {
                title: "Add combine"
                iconSource: "qrc:///icons/join_64.png"
                MenuItem { action: addMinimum }
                MenuItem { action: addMaximum }
                MenuItem { action: addInhibition }
            }

            Menu {
                title: "Add transform"
                iconSource: "qrc:///icons/trafo_64.png"
                MenuItem { action: addAffine }
            }

            Menu {
                title: "Add output"
                iconSource: "qrc:///icons/output_64.png"
                MenuItem { action: addMarginalPlot }
                MenuItem { action: addScatterPlot }
                MenuItem { action: addKDEPlot }
            }
        }
        Menu {
            title: "View"
            MenuItem { action: zoomIn }
            MenuItem { action: zoomOut }
        }
        Menu {
            title: "Help"
            MenuItem { action: showAbout }
            MenuItem { action: showManual }
        }
    }

    toolBar: ToolBar {
        RowLayout {
            ToolButton { action: verifyNet }
            ToolButton { action: runNet }
            //ToolBarSeparator {}
            ToolButton {
                text: "Add Variable"
                iconSource: "qrc:///icons/var_64.png"
                menu: Menu {
                    MenuItem { action: addStimulus }
                    MenuSeparator {}
                    MenuItem { action: addConstant }
                    MenuItem { action: addUniformVariable }
                    MenuItem { action: addNormalVariable }
                    MenuItem { action: addGammaVariable }
                    MenuItem { action: addInverseGammaVariable }
                    MenuItem { action: addWeibullVariable }
                    MenuSeparator {}
                    MenuItem { action: addCompoundNormalVariable }
                    MenuItem { action: addCompoundGammaVariable }
                    MenuItem { action: addCompoundInverseGammaVariable }
                    MenuItem { action: addCompoundWeibullVariable }
                }
            }
            ToolButton {
                text: "Add Process"
                iconSource: "qrc:///icons/proc_64.png"
                menu: Menu {
                    MenuItem { action: addFixedDelay }
                    MenuItem { action: addGammaProcess }
                    MenuItem { action: addInverseGammaProcess }
                    MenuItem { action: addWeibullProcess }
                    MenuSeparator {}
                    MenuItem { action: addRandomDelay }
                    MenuItem { action: addCompoundGammaProcess }
                    MenuItem { action: addCompoundInverseGammaProcess }
                    MenuItem { action: addCompoundWeibullProcess }
                }
            }
            ToolButton {
                text: "Add process"
                iconSource: "qrc:///icons/proc_64.png"
                menu: Menu {
                    MenuItem { action: addFixedDelay }
                    MenuItem { action: addGammaProcess }
                    MenuItem { action: addInverseGammaProcess }
                    MenuItem { action: addWeibullProcess }
                    MenuSeparator {}
                    MenuItem { action: addRandomDelay }
                    MenuItem { action: addCompoundGammaProcess }
                    MenuItem { action: addCompoundInverseGammaProcess }
                    MenuItem { action: addCompoundWeibullProcess }
                }
            }
            ToolButton {
                text: "Add combine"
                iconSource: "qrc:///icons/join_64.png"
                menu: Menu {
                    MenuItem { action: addMinimum }
                    MenuItem { action: addMaximum }
                    MenuItem { action: addInhibition }
                }
            }
            ToolButton {
                text: "Add transform"
                iconSource: "qrc:///icons/trafo_64.png"
                menu: Menu {
                    MenuItem { action: addAffine }
                }
            }
            ToolButton {
                text: "Add output"
                iconSource: "qrc:///icons/output_64.png"
                menu: Menu {
                    MenuItem { action: addMarginalPlot }
                    MenuItem { action: addScatterPlot }
                    MenuItem { action: addKDEPlot }
                }
            }
            //ToolBarSeparator {}
            ToolButton { action: zoomIn }
            ToolButton { action: zoomOut }
            //ToolBarSeparator {}
            ToolButton { action: deleteItem }
        }
    }

    Network {
        id: networkEdit
    }
}
