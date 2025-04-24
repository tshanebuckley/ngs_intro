#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path "${greeting}-output.txt"

    script:
    """
    echo '$greeting' > '$greeting-output.txt'
    """
}

process convertToUpper {

    publishDir 'results', mode: 'copy'

    input:
        path input_file

    output:
        path "UPPER-${input_file}.txt"

    script:
    """
    cat $input_file | tr '[a-z]' '[A-Z]' > UPPER-${input_file}.txt
    """
}

process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
        path input_files
        val batch_name

    output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        val count_greetings, emit: count

    script:
    count_greetings = input_files.size()
    """
    cat $input_files > COLLECTED-${batch_name}-output.txt
    """
}


/*
 * Pipeline parameters
 */
params.greeting = 'greetings.csv'
params.batch_name = 'test'

workflow {

    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath(params.greeting)
                        .splitCsv()
                        .map { line -> line[0] }

    // emit a greeting
    sayHello(greeting_ch)

    convertToUpper(sayHello.out)

    collectGreetings(convertToUpper.out.collect(), params.batch_name) 

    collectGreetings.out.count.view { count -> "$count greeting(s)"}
}
