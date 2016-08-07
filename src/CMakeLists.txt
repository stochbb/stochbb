SET(stochbb_SOURCES main.cc assembler.cc
    qcustomplot.cc qnetview.cc mainwindow.cc neteditwidget.cc nodes.cc edge.cc network.cc
    plotwindow.cc logwindow.cc parameter.cc)
SET(stochbb_MOC_HEADERS
    qcustomplot.hh qnetview.hh mainwindow.hh neteditwidget.hh nodes.hh edge.hh network.hh
    plotwindow.hh logwindow.hh parameter.hh)
SET(stochbb_HEADERS assembler.hh ${stochbb_MOC_HEADERS})

SET(LANGUAGES de_DE)
SET(LANGUAGE_TS_FILES)
FOREACH(LANGUAGE ${LANGUAGES})
  SET(TS_FILE "${CMAKE_SOURCE_DIR}/shared/i18n/stochbb_${LANGUAGE}.ts")
  SET(LANGUAGE_TS_FILES ${LANGUAGE_TS_FILES} ${TS_FILE})
  SET_SOURCE_FILES_PROPERTIES(${TS_FILE} PROPERTIES OUTPUT_LOCATION
    "${CMAKE_SOURCE_DIR}/shared/i18n")
ENDFOREACH(LANGUAGE ${LANGUAGES})
IF(LANGUAGE_TS_FILES)
  QT5_CREATE_TRANSLATION(stochbb_QM_FILES ${stochbb_SOURCES} ${LANGUAGE_TS_FILES})
ENDIF(LANGUAGE_TS_FILES)

qt5_wrap_cpp(stochbb_MOC_SOURCES ${stochbb_MOC_HEADERS})
qt5_add_resources(stochbb_RCC_SOURCES ../shared/resources.qrc)

add_executable(stochbb ${stochbb_SOURCES} ${stochbb_MOC_SOURCES} ${stochbb_RCC_SOURCES} ${stochbb_QM_FILES})
target_link_libraries(stochbb ${LIBS})

install(TARGETS stochbb DESTINATION "bin/")

if (UNIX AND NOT APPLE)
 install(FILES ../shared/stochbb.desktop DESTINATION share/applications)
 install(FILES ../shared/stochbb-icon.png DESTINATION share/icons)
endif (UNIX AND NOT APPLE)
