version 1.0

# WORKFLOW DEFINITION 
workflow TestWF {
  input {
    String sample_basename
  }
  
  call FileToArray
  
  scatter (chromosome in FileToArray.scatter_regions ) {
    call DoSomethingWithChromosome {
      input:
        chromosome=chromosome,
        sample_basename=sample_basename
    }

  }

}

task FileToArray {
    String fname = "Test_chromosome_interval_list.txt"
    String github_fname = "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/softsearch/Test_chromosome_interval_list.txt"
    command <<<
        # Get chromosome interval list file
        ## wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/softsearch/Test_chromosome_interval_list.txt
        wget  ~{github_fname}
        # Re-attempt the wget without the https_proxy set
        unset https_proxy 
        ##wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/softsearch/Test_chromosome_interval_list.txt
        wget  ~{github_fname}     
        ## cat Test_chromosome_interval_list.txt
        cat ~{fname}
    >>>
      
    runtime {
        docker: "alesmaver/softsearch"
        runtime_minutes: 5
    }
    output {
        Array[String] scatter_regions = read_lines(stdout())
    }
}


task DoSomethingWithChromosome {
  input {
    String chromosome
    String sample_basename
  }

  command <<<
    echo sample_basename: ~{sample_basename}
    echo chromosome: ~{chromosome}
  >>>

  runtime {
    docker: "alesmaver/softsearch"
    runtime_minutes: 5
  }
  output {
    #File outfile = stdout()
  }
}

