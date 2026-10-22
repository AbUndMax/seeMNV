package com.github.abundmax.seemnv.core.fieldArray;

import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Properties;

public record FieldConfig(String id, String number, FieldType fieldType, FieldAggregationStrategy aggregationStrategy) {

    private static final Path AGGREGATION_PROPERTIES_PATH = Paths.get("config", "aggregation.properties");
    private static final Properties CONFIG = new Properties();

    static {
        try (Reader reader = Files.newBufferedReader(AGGREGATION_PROPERTIES_PATH)) {

            CONFIG.load(reader);

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static FieldConfig buildFromHeaderString(String headerLineFormat) {
        // eg.: ##FORMAT=<ID=AD,Number=R,Type=Integer>
        headerLineFormat = headerLineFormat.split("=<")[1];
        // ["##FORMAT", "ID=AD,Number=R,Type=Integer>"]
        String[] split = headerLineFormat.split(",");
        // ["ID=AD", "Number=R", "Type=Integer>"]

        String id = split[0].split("=")[1]; // "AD"

        if (id.equalsIgnoreCase("GT")) {
            return null;
        }

        String number = split[1].split("=")[1]; // "R"
        String type = split[2].split("=")[1]; // "Integer"
        FieldType fieldType = FieldType.fromString(type);

        String definedStrategy = CONFIG.getProperty(id);
        if (definedStrategy == null) {
            throw new IllegalArgumentException("No valid FORMAT ID specified in aggregation.properties for ID: " + id);
        }

        FieldAggregationStrategy aggregationStrategy = FieldAggregationStrategy.fromString(definedStrategy);

        return new FieldConfig(id, number, fieldType, aggregationStrategy);
    }
}