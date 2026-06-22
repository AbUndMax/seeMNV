package com.github.abundmax.seemnv.cli;


public class Main {

    public static void main(String[] args) {
        System.out.println("Hello from seeMNV!");

        // args enthält die Kommandozeilenargumente
        for (String arg : args) {
            System.out.println(arg);
        }
    }
}