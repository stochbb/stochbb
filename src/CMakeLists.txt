SET(stochbb_SOURCES main.cc assembler.cc
    qcustomplot.cc qnetview.cc mainwindow.cc neteditwidget.cc nodes.cc edge.cc network.cc
    plotwindow.cc logwindow.cc parameter.cc)
SET(stochbb_MOC_HEADERS
    qcustomplot.hh qnetview.hh mainwindow.hh neteditwidget.hh nodes.hh edge.hh network.hh
    plotwindow.hh logwindow.hh parameter.hh)
SET(stochbb_HEADERS assembler.hh ${stochbb_MOC_HEADERS})

qt5_wrap_cpp(stochbb_MOC_SOURCES ${stochbb_MOC_HEADERS})
qt5_add_resources(stochbb_RCC_SOURCES ../shared/resources.qrc)

add_executable(stochbb ${stochbb_SOURCES} ${stochbb_MOC_SOURCES} ${stochbb_RCC_SOURCES})
target_link_libraries(stochbb ${LIBS} libsbb)

install(TARGETS stochbb DESTINATION "bin/")
