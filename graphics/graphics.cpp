#include <QApplication>
#include <QMainWindow>
#include <QChartView>
#include <QScatterSeries>
#include <QValueAxis>
#include <QLineSeries>
#include <QFile>
#include <boost/json.hpp>
#include <vector>
#include <iostream>
#include "graphics.h"

using namespace QtCharts;
using namespace boost::json;

// Function to visualize points and constraints in a Qt window
void visualizePoints(const std::vector<std::pair<double, double>>& points, const std::vector<std::pair<int, int>>& constraints) {
    // Create a main window
    QMainWindow mainWindow;

    // Create a scatter series to hold the points
    QScatterSeries* series = new QScatterSeries();
    series->setMarkerSize(10); // Adjust the marker size for better visibility
    for (const auto& point : points) {
        series->append(point.first, point.second);
    }

    // Create a line series for the constraints
    QLineSeries* lineSeries = new QLineSeries();
    
    // Draw lines for each constraint
    for (const auto& constraint : constraints) {
        int point1_idx = constraint.first;
        int point2_idx = constraint.second;

        // Ensure the indices are valid
        if (point1_idx < points.size() && point2_idx < points.size()) {
            lineSeries->append(points[point1_idx].first, points[point1_idx].second);
            lineSeries->append(points[point2_idx].first, points[point2_idx].second);
            lineSeries->append(QPointF(points[point1_idx].first, points[point1_idx].second)); // Close the line segment
            lineSeries->append(QPointF(points[point2_idx].first, points[point2_idx].second)); // Close the line segment
        }
    }

    // Create a chart and add the series
    QChart* chart = new QChart();
    chart->addSeries(series);
    chart->addSeries(lineSeries); // Add the constraints as a separate series
    chart->setTitle("Point and Constraint Visualization");
    chart->createDefaultAxes();

    // Set the axes titles
    QValueAxis* axisX = new QValueAxis();
    axisX->setTitleText("X Axis");
    axisX->setRange(0, 10000); // Adjust the range as needed

    QValueAxis* axisY = new QValueAxis();
    axisY->setTitleText("Y Axis");
    axisY->setRange(0, 10000); // Adjust the range as needed

    chart->addAxis(axisX, Qt::AlignBottom);
    chart->addAxis(axisY, Qt::AlignLeft);
    series->attachAxis(axisX);
    series->attachAxis(axisY);
    lineSeries->attachAxis(axisX);
    lineSeries->attachAxis(axisY);

    // Create a chart view to display the chart
    QChartView* chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);
    mainWindow.setCentralWidget(chartView);
    mainWindow.resize(800, 600);
    mainWindow.show();

    // Start the Qt application event loop
    QApplication::exec();
}



