#!/usr/bin/env python
# coding: utf-8

# In[1]:


"Basic Quantum Network Example: components, nodes, and connections"
#####################################################################
import netsquid as ns

node_alice=ns.nodes.Node(name='Alice')
node_bob=ns.nodes.Node(name='Bob')
#create Alice and Bob nodes

#Now we create a 'delay model' of our channel using the avg speed (0.75c) at which information travels through the channel, plus standard dev. (0.05c)
class alice_bob_delay_model(ns.components.models.DelayModel):
    def __init__(self, speed_of_light_fraction=0.75, standard_deviation=0.05):
        super().__init__() 
        #I think this just ensures that my_delay_model inherits characteristics from netsquid's DelayModel class
        self.properties["speed"]=speed_of_light_fraction *3e5 #in km/s
        self.properties["std"]=standard_deviation
        self.required_properties=['length'] #in km
    def generate_delay(self,**kwargs):
        avg_speed=self.properties["speed"]
        std=self.properties["std"]
        #Below is directly from the netsquid docs: 'rng' is a property containing random number gen.
        #We use this to randomly select a speed from a normal distribution
        speed=self.properties["rng"].normal(avg_speed,avg_speed*std)
        delay=1e9*kwargs['length'] /speed #delay is in nanoseconds, hence the 1e9
        return delay

#Now create our quantum channel: 

distance=2.74/1000 #We need to set the length of the channel (2.74 m) in km
my_delay_model=alice_bob_delay_model() #establish the current delay model to use in the alice/bob channel

channel_1=ns.components.QuantumChannel(name='qchannel[alice to bob]',
                                      length=distance,
                                      models={"delay_model":my_delay_model})
channel_2=ns.components.QuantumChannel(name='qchannel[bob to alice]',
                                      length=distance,
                                      models={"delay_model":my_delay_model})
#Channel 1: Alice to Bob
#Channel 2: Bob to Alice

#Now we create a 'connection' to wrap the channels into a component which connects alice and bob's nodes
connection=ns.nodes.DirectConnection(name='conn[alice|bob]',
                                    channel_AtoB=channel_1,
                                    channel_BtoA=channel_2)
#connect the nodes:
node_alice.connect_to(remote_node=node_bob, 
                      connection=connection, 
                      local_port_name='qubitIO',
                      remote_port_name='qubitIO')


# In[2]:


"Protocols: Bob waits for Alice to send a qubit. Upon recieving a qubit, Bob measures it in some basis and then sends it back."
#NOTE: for this cell to run properly, clear all variables before re-running

class alice_bob_protocol(ns.protocols.NodeProtocol):
    def __init__(self, node, observable, qubit=None):
        super().__init__(node)
        self.observable=observable
        self.qubit=qubit
        self.basis=["|0>","|1>"] if observable==ns.Z else ["|+>","|->"] #this is just for printing I guess
        
    def run(self):
        if self.qubit is not None:
            self.node.ports["qubitIO"].tx_output(self.qubit) #??????
        while True:
            yield self.await_port_input(self.node.ports["qubitIO"]) #wait until qubit arrives at the port
            ###############################
            #I don't understand this section.
            message=self.node.ports["qubitIO"].rx_input()
            qubit=message.items[0]
            meas,prob=ns.qubits.measure(qubit,observable=self.observable)
            print(f"{ns.sim_time():5.1f}: {self.node.name} measured "
                  f"{self.basis[meas]} with probability {prob:.2f}")
            self.node.ports["qubitIO"].tx_output(qubit)
            #######################################

#assign the protocol to alice and bob's nodes
qubits=ns.qubits.create_qubits(1)
alice_protocol=alice_bob_protocol(node_alice, observable=ns.Z, qubit=qubits[0])
bob_protocol=alice_bob_protocol(node_bob, observable=ns.X)

#run the simulation
alice_protocol.start()
bob_protocol.start()
run_stats=ns.sim_run(duration=100)
print(run_stats)

#you should see Alice measuring |0> or |1> with 50/50 probability and Bob measuring |+> or |-> with 50/50 probability.


# In[3]:


print("{}".format(ns.examples.teleportation.__file__))


# In[ ]:




