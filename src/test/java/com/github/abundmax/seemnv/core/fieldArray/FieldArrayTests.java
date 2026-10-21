package com.github.abundmax.seemnv.core.fieldArray;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;

public class FieldArrayTests {

    @Test
    public void testIntegerFieldArrayCreationWithSingleValues() {
        FieldConfig fieldConfig = createConfig("GQ", "1", "Integer");

        FieldArray<?> fieldArray = FieldArrayFabricator.get(fieldConfig)
                .apply(new String[]{"1", "3"});

        assertSame(fieldConfig, fieldArray.fieldConfig);
        assertEquals(2, fieldArray.length);
        assertArrayEquals(new Integer[]{1}, fieldArray.values[0]);
        assertArrayEquals(new Integer[]{3}, fieldArray.values[1]);
    }

    @Test
    public void testIntegerFieldArrayCreationWithMultipleValuesPerSample() {
        FieldConfig fieldConfig = createConfig("GQ", "2", "Integer");

        FieldArray<?> fieldArray = FieldArrayFabricator.get(fieldConfig)
                .apply(new String[]{"1,2", "3,4"});

        assertSame(fieldConfig, fieldArray.fieldConfig);
        assertEquals(2, fieldArray.length);
        assertArrayEquals(
                new Integer[][]{
                        new Integer[]{1, 2},
                        new Integer[]{3, 4}
                },
                fieldArray.values
        );
    }

    @Test
    public void testFloatFieldArrayCreation() {
        FieldConfig fieldConfig = createConfig("AF", "A", "Float");

        FieldArray<?> fieldArray = FieldArrayFabricator.get(fieldConfig)
                .apply(new String[]{"0.1", "0.25", "0.004831"});

        assertSame(fieldConfig, fieldArray.fieldConfig);
        assertEquals(3, fieldArray.length);
        assertArrayEquals(new Float[]{0.1f}, fieldArray.values[0]);
        assertArrayEquals(new Float[]{0.25f}, fieldArray.values[1]);
        assertArrayEquals(new Float[]{0.004831f}, fieldArray.values[2]);
    }

    @Test
    public void testFloatFieldArrayCreationWithMultipleValuesPerSample() {
        FieldConfig fieldConfig = createConfig("AF", "A", "Float");

        FieldArray<?> fieldArray = FieldArrayFabricator.get(fieldConfig)
                .apply(new String[]{"0.1,0.2", "0.3,0.4"});

        assertSame(fieldConfig, fieldArray.fieldConfig);
        assertEquals(2, fieldArray.length);
        assertArrayEquals(new Float[]{0.1f, 0.2f}, fieldArray.values[0]);
        assertArrayEquals(new Float[]{0.3f, 0.4f}, fieldArray.values[1]);
    }

    @Test
    public void testBooleanFieldArrayCreation() {
        FieldConfig fieldConfig = createConfig("GQ", "1", "Boolean");

        FieldArray<?> fieldArray = FieldArrayFabricator.get(fieldConfig)
                .apply(new String[]{"true", "false", "true"});

        assertSame(fieldConfig, fieldArray.fieldConfig);
        assertEquals(3, fieldArray.length);
        assertArrayEquals(new Boolean[]{true}, fieldArray.values[0]);
        assertArrayEquals(new Boolean[]{false}, fieldArray.values[1]);
        assertArrayEquals(new Boolean[]{true}, fieldArray.values[2]);
    }

    @Test
    public void testStringFieldArrayCreation() {
        FieldConfig fieldConfig = createConfig("GQ", "1", "String");

        FieldArray<?> fieldArray = FieldArrayFabricator.get(fieldConfig)
                .apply(new String[]{"AA", "BB", "CC"});

        assertSame(fieldConfig, fieldArray.fieldConfig);
        assertEquals(3, fieldArray.length);
        assertArrayEquals(new String[]{"AA"}, fieldArray.values[0]);
        assertArrayEquals(new String[]{"BB"}, fieldArray.values[1]);
        assertArrayEquals(new String[]{"CC"}, fieldArray.values[2]);
    }

    @Test
    public void testStringFieldArrayCreationWithMultipleValuesPerSample() {
        FieldConfig fieldConfig = createConfig("GQ", "2", "String");

        FieldArray<?> fieldArray = FieldArrayFabricator.get(fieldConfig)
                .apply(new String[]{"AA,BB", "CC,DD"});

        assertSame(fieldConfig, fieldArray.fieldConfig);
        assertEquals(2, fieldArray.length);
        assertArrayEquals(new String[]{"AA", "BB"}, fieldArray.values[0]);
        assertArrayEquals(new String[]{"CC", "DD"}, fieldArray.values[1]);
    }

    @Test
    public void testCombinationIntUsesElementWiseMinimum() {
        FieldConfig fieldConfig = createConfig("GQ", "2", "Integer");

        FieldArray<?> fieldArray = FieldArrayFabricator.get(fieldConfig)
                .apply(new String[]{"1,2", "3,4"});

        FieldArray<?> fieldArray2 = FieldArrayFabricator.get(fieldConfig)
                .apply(new String[]{"5,1", "2,8"});

        FieldArray<?> fieldArray3 = FieldArrayFabricator.get(fieldConfig)
                .apply(new String[]{"0,5", "1,5"});

        FieldArray<?> combined = FieldArray.combineWildcard(
                fieldArray,
                fieldArray2,
                fieldArray3
        );

        assertSame(fieldConfig, combined.fieldConfig);
        assertEquals(2, combined.length);
        assertArrayEquals(new Integer[]{0, 1}, combined.values[0]);
        assertArrayEquals(new Integer[]{1, 4}, combined.values[1]);
    }

    @Test
    public void testCombinationFloatUsesElementWiseMinimum() {
        FieldConfig fieldConfig = createConfig("GQ", "2", "Float");

        FieldArray<?> fieldArray = FieldArrayFabricator.get(fieldConfig)
                .apply(new String[]{"1.3,2.5", "3.6,4.9"});

        FieldArray<?> fieldArray2 = FieldArrayFabricator.get(fieldConfig)
                .apply(new String[]{"5.7,8.7", "2.9,8.8"});

        FieldArray<?> fieldArray3 = FieldArrayFabricator.get(fieldConfig)
                .apply(new String[]{"7.7,5.9", "9.7,5.4"});

        FieldArray<?> combined = FieldArray.combineWildcard(
                fieldArray,
                fieldArray2,
                fieldArray3
        );

        assertSame(fieldConfig, combined.fieldConfig);
        assertEquals(2, combined.length);
        assertArrayEquals(new Float[]{1.3f, 2.5f}, combined.values[0]);
        assertArrayEquals(new Float[]{2.9f, 4.9f}, combined.values[1]);
    }

    @Test
    public void testCombinationUsesElementWiseMaximumForPL() {
        FieldConfig fieldConfig = createConfig("PL", "G", "Integer");

        FieldArray<?> fieldArray = FieldArrayFabricator.get(fieldConfig)
                .apply(new String[]{"0,100", "10,20"});

        FieldArray<?> fieldArray2 = FieldArrayFabricator.get(fieldConfig)
                .apply(new String[]{"5,80", "30,15"});

        FieldArray<?> fieldArray3 = FieldArrayFabricator.get(fieldConfig)
                .apply(new String[]{"2,120", "25,40"});

        FieldArray<?> combined = FieldArray.combineWildcard(
                fieldArray,
                fieldArray2,
                fieldArray3
        );

        assertSame(fieldConfig, combined.fieldConfig);
        assertEquals(2, combined.length);
        assertArrayEquals(new Integer[]{5, 120}, combined.values[0]);
        assertArrayEquals(new Integer[]{30, 40}, combined.values[1]);
    }

    @Test
    public void testCombinationWithSingleFieldArrayReturnsEquivalentValues() {
        FieldConfig fieldConfig = createConfig("GQ", "2", "Integer");

        FieldArray<?> fieldArray = FieldArrayFabricator.get(fieldConfig)
                .apply(new String[]{"1,2", "3,4"});

        FieldArray<?> combined = FieldArray.combineWildcard(fieldArray);

        assertSame(fieldConfig, combined.fieldConfig);
        assertEquals(2, combined.length);
        assertArrayEquals(new Integer[]{1, 2}, combined.values[0]);
        assertArrayEquals(new Integer[]{3, 4}, combined.values[1]);
    }

    @Test
    public void testCombineWildcardThrowsForNoFieldArrays() {
        assertThrows(
                IllegalArgumentException.class,
                FieldArray::combineWildcard
        );
    }

    @Test
    public void testCombineWildcardThrowsForNullFieldArrays() {
        assertThrows(
                IllegalArgumentException.class,
                () -> FieldArray.combineWildcard((FieldArray<?>[]) null)
        );
    }

    @Test
    public void testCombineWildcardThrowsForDifferentSampleLengths() {
        FieldConfig fieldConfig = createConfig("GQ", "1", "Integer");

        FieldArray<?> fieldArray = FieldArrayFabricator.get(fieldConfig)
                .apply(new String[]{"1", "2"});

        FieldArray<?> fieldArray2 = FieldArrayFabricator.get(fieldConfig)
                .apply(new String[]{"1", "2", "3"});

        assertThrows(
                IllegalStateException.class,
                () -> FieldArray.combineWildcard(fieldArray, fieldArray2)
        );
    }

    @Test
    public void testCombineWildcardThrowsForDifferentFieldConfigInstances() {
        FieldConfig fieldConfig = createConfig("GQ", "1", "Integer");
        FieldConfig otherFieldConfig = createConfig("GQ", "1", "Integer");

        FieldArray<?> fieldArray = FieldArrayFabricator.get(fieldConfig)
                .apply(new String[]{"1", "2"});

        FieldArray<?> fieldArray2 = FieldArrayFabricator.get(otherFieldConfig)
                .apply(new String[]{"3", "4"});

        assertThrows(
                IllegalArgumentException.class,
                () -> FieldArray.combineWildcard(fieldArray, fieldArray2)
        );
    }

    @Test
    public void testCombineWildcardThrowsForDifferentValueTypes() {
        FieldConfig intFieldConfig = createConfig("GQ", "1", "Integer");
        FieldConfig floatFieldConfig = createConfig("AF", "1", "Float");

        FieldArray<?> intFieldArray = FieldArrayFabricator.get(intFieldConfig)
                .apply(new String[]{"1", "2"});

        FieldArray<?> floatFieldArray = FieldArrayFabricator.get(floatFieldConfig)
                .apply(new String[]{"1.5", "2.5"});

        assertThrows(
                IllegalArgumentException.class,
                () -> FieldArray.combineWildcard(intFieldArray, floatFieldArray)
        );
    }

    private FieldConfig createConfig(String id, String number, String type) {
        return new FieldConfig(
                "##FORMAT=<ID=" + id +
                        ",Number=" + number +
                        ",Type=" + type +
                        ",Description=\"Test field\">"
        );
    }
}