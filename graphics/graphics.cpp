#include <QApplication>
#include <QMainWindow>
#include <QChartView>
#include <QScatterSeries>
#include <QValueAxis>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QFile>
#include <boost/json.hpp>
#include <vector>
#include <iostream>
#include "graphics.h"

using namespace QtCharts;
using namespace boost::json; // Use Boost JSON namespace

// Function to read points from JSON file
std::vector<std::pair<double, double>> readPointsFromJson(const std::string& filename) {
    std::vector<std::pair<double, double>> points;

    // Open the JSON file
    QFile file(QString::fromStdString(filename));
    if (!file.open(QIODevice::ReadOnly)) {
        std::cerr << "Could not open the file: " << filename << std::endl;
        return points; // Return empty vector on error
    }

    // Read the file into a JSON document
    QByteArray jsonData = file.readAll();
    file.close();

    // Parse JSON using Boost.JSON
    try {
        value json_value = parse(std::string(jsonData.constData()));
        const auto& points_array = json_value.as_object()["points"].as_array();

        // Read points from JSON
        for (const auto& item : points_array) {
            const auto& obj = item.as_object();
            double x = obj.at("x").as_double();
            double y = obj.at("y").as_double();
            points.emplace_back(x, y);
        }
    } catch (const std::exception& e) {
        std::cerr << "Error parsing JSON: " << e.what() << std::endl;
    }

    return points;
}

// Function to visualize points in a Qt window
void visualizePoints(const std::vector<std::pair<double, double>>& points) {
    // Create a main window
    QMainWindow mainWindow;

    // Create a scatter series to hold the points
    QScatterSeries* series = new QScatterSeries();
    for (const auto& point : points) {
        series->append(point.first, point.second);
    }

    // Create a chart and add the series
    QChart* chart = new QChart();
    chart->addSeries(series);
    chart->setTitle("Point Visualization");
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


    // Create a chart view to display the chart
    QChartView* chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);
    mainWindow.setCentralWidget(chartView);
    mainWindow.resize(800, 600);
    mainWindow.show();

    // Start the Qt application event loop
    QApplication::exec();
}




