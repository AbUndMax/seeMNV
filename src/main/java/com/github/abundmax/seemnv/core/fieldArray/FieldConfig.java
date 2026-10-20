package com.github.abundmax.seemnv.core.fieldArray;

import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Properties;

public class FieldConfig {

    private static final Path AGGREGATION_PROPERTIES_PATH = Paths.get("config", "aggregation.properties");
    private static final Properties CONFIG = new Properties();

    static {
        try (Reader reader = Files.newBufferedReader(AGGREGATION_PROPERTIES_PATH)) {

            CONFIG.load(reader);

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public final String id;
    public final String number;
    public final FieldType fieldType;
    public final FieldAggregationStrategy aggregationStrategy;

    public FieldConfig(String headerLineFormat) {
        // eg.: ##FORMAT=<ID=AD,Number=R,Type=Integer>
        headerLineFormat = headerLineFormat.split("=<")[1];
        // ["##FORMAT", "ID=AD,Number=R,Type=Integer>"]
        String[] split = headerLineFormat.split(",");
        // ["ID=AD", "Number=R", "Type=Integer>"]

        id = split[0].split("=")[1]; // "AD"
        number = split[1].split("=")[1]; // "R"
        String type = split[2].split("=")[1]; // "Integer"
        fieldType = FieldType.fromString(type);

        String definedStrategy = CONFIG.getProperty(id);
        if (definedStrategy == null) {
            throw new IllegalArgumentException("No valid FORMAT ID specified in aggregation.properties for ID: " + id);
        }

        aggregationStrategy = FieldAggregationStrategy.fromString(definedStrategy);
    }
}
