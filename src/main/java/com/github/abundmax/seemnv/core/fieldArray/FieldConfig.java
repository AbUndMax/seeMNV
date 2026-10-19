package com.github.abundmax.seemnv.core.fieldArray;

import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;
import java.util.Properties;

public class FieldConfig {

    private static final Map<String, Class<?>> STRING_TO_CLASS = Map.of(
            "string", String.class,
            "integer", Integer.class,
            "float", Float.class,
            "boolean", Boolean.class
    );

    private static final Path AGGREGATION_PROPERTIES_PATH = Paths.get("config", "aggregation.properties");
    private static final Properties CONFIG = new Properties();

    static {
        try (Reader reader = Files.newBufferedReader(AGGREGATION_PROPERTIES_PATH)) {

            CONFIG.load(reader);

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public final String ID;
    public final String NUMBER;
    public final FieldType FIELD_TYPE;
    public final FieldAggregationStrategy AGGREGATION_STRATEGY;

    public FieldConfig(String headerLineFormat) {
        // eg.: ##FORMAT=<ID=AD,Number=R,Type=Integer>
        headerLineFormat = headerLineFormat.split("=<")[1];
        // ["##FORMAT", "ID=AD,Number=R,Type=Integer>"]
        String[] split = headerLineFormat.split(",");
        // ["ID=AD", "Number=R", "Type=Integer>"]

        ID = split[0].split("=")[1]; // "AD"
        NUMBER = split[1].split("=")[1]; // "R"
        String type = split[2].split("=")[1]; // "Integer"
        FIELD_TYPE = FieldType.fromString(type);

        String definedStrategy = CONFIG.getProperty(ID);
        if (definedStrategy == null) {
            throw new IllegalArgumentException("No valid FORMAT ID specified in aggregation.properties for ID: " + ID);
        }

        AGGREGATION_STRATEGY = FieldAggregationStrategy.fromString(definedStrategy);
    }
}
